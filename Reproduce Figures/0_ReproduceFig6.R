#################################################################################################
##### PaintSims_Study.R
##### 26 November 2024
##### Simulations to Run on a Cluster
#################################################################################################

### Works best when run on a cluster
### This code reproduces the frame in Figure 6 associated with n.sources defined below

#################################################################################################
##### Initialize Workspace 
#################################################################################################
### Start from a clean workspace 
rm(list=ls())
library(dplyr) 
library(ggplot2)

#################################################################################################
###                                       Info to Change                                        #
#################################################################################################
### Define paths                                                                                #
# Set working directory - this is the folder where your data files are saved                    #
setwd('/Users/mausdemore/Documents/GitHub/kenary/Technometrics_Submit')                         #                                       
                                                                                                #
### Define the number of cores that should be used                                              #
n.cores = parallel::detectCores()-2                                                             #
                                                                                                #
### Define the number of paint source combinations to consider                                  # 
# n.combos= 100 for full set of simulations considered in paper, can consider n.combos<100 for  #
# a quicker study                                                                               #
n.combos = 25                                                                                   #
                                                                                                #
### Define the parameters for the study                                                         #
# Define the number of sources to be studied - n = 2,3,4,5,6,7, or 8 in paper                   #
n.sources = 4                                                                                   #
# Define the number of control objects - n.control.samps=c(3:9) considered in paper             #
n.control.samps = c(3:9)                                                                        #
# Define the number of test objects - n.test.samps=3 for paper                                  #
n.test.samps = 3                                                                                #
#################################################################################################
###                       !!! Do not change anything below this line !!!                        #
#################################################################################################


#################################################################################################
##### Load the Data 
#################################################################################################
### Source Function File 
source("KeNary_Funs.R")                                                                         

### Load paint data 
load('Sim_Spectra.rda')
# Get wave numbers from data (x-axis)
wav.num = sim.spectra[[1]][,1]

### Load source information
source.grid = as.matrix(fread('source_grid.csv'), nrow=100, ncol=10)
source.grid = source.grid[1:n.combos,]
#################################################################################################


#################################################################################################
##### Set up Model 
#################################################################################################
### Define relevant spectral info
# How many basis functions do we want to use to define the spectra?
p.basis = 300
# How many points do we want to use to evaluate the pseudo-spectra?
p.eval = 1000
# At which specific points do we want to evaluate the spectra?
xs.eval = wav.num[round(seq(1,length(wav.num),length=p.eval))]
# What lag do we want to consider for the cross-correlation function?
max.lag = 10

### Define kernel function for spectra
spectkern = spectral.kernel(xs.eval, max.lag, wav.num)

### Define list to store data 
FTIR.sims.results = list()
FTIR.sims.results[[1]]=FTIR.sims.results[[2]]=list()
for(cs in n.control.samps){
     ### Status update
     print(cs)
     s = n.sources
     
     ### Update n.control.samps if not a vector 
     cs <-unique(cs)
     if(length(cs)==1){
          cs = rep(cs, s)
     }
     
     ### Define object pairwise comparisons for parameter estimation
     # Get an initial grid of pairwise comparisons 
     pairs.all = t(combn(1:(sum(cs)+n.test.samps),2))
     # Define indices to distinguish between trace and control objects - which.prop here doesn't matter 
     ix.all = get.indices(n.sources=s, n.sizes=cs, 
                           n.unknown=n.test.samps, which.prop=paste("H", s, sep=''))
     # Grab control pairs 
     control.ix = which(unlist(ix.all[1:(s+choose(s, 2))])%in% ix.all$control.IX)
     control.pairs = pairs.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], ]
     
     
     ### Create P matrices
     P.mat.all = make.P(sum(cs)+n.test.samps)
     P.mat.control = P.mat.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], 
                                1:sum(cs)] 
     PP.t.control = P.mat.control %*% t(P.mat.control)
     
     
     ### Prepare Adaptive MH sampler
     # Define starting point for sampler 
     theta0 = c(rep(50,choose(s,2)+s), rep(5, choose(s,2)+s), 0.1) 
     
     ### Get performance of Adaptive MH versus SVM for paint data
     registerDoParallel(cores=n.cores)
     FTIR.sims.results[[unique(cs)]] = foreach(i=1:nrow(source.grid)) %dopar%
                                FTIR.sims(sims.info=source.grid[i,], dat.ls=sim.spectra, n.sources=s, n.control=cs,
                                          n.test=n.test.samps, control.pairs=control.pairs, 
                                          PP.t.control=PP.t.control, control.ix=control.ix, P.mat.all=P.mat.all, 
                                          theta0=theta0, package.vect=package.vect, kernfun=spectkern, wav.num=wav.num, 
                                          SVM.analysis=TRUE, seed=0806, n.samps=5e05)
     stopImplicitCluster()
}

### Recreate Part of Figure 6
# Get SVM versus KeNary Performance
n.source = s
KBSC.rate = SVM.rate = list()
for(i in 3:length(FTIR.sims.results)){
        SVM.rate[[i]] = matrix(0, nrow=length(FTIR.sims.results[[i]]), ncol=3)
        KBSC.rate[[i]] = sum(unlist(lapply(FTIR.sims.results[[i]], function(x) which.max(x$BF)==x$n.sources)))/length(FTIR.sims.results[[i]])
        for(j in 1:length(FTIR.sims.results[[i]])){
                tmp.col = apply(FTIR.sims.results[[i]][[j]]$p.SVM, 1, which.max)
                SVM.rate[[i]][j,] = names(sapply(1:3, function(x, sims.results, tmp.col) 
                        sims.results[x, tmp.col[x]], FTIR.sims.results[[i]][[j]]$p.SVM, tmp.col)) == source.grid[j, n.source]
        }
}

### Look at Kenary Success Rate for 2 classes as number of control objects varies 3:10
KBSC.df = do.call('rbind', lapply(KBSC.rate, unlist))
rownames(KBSC.df) = paste(n.control.samps, 'control obs')
colnames(KBSC.df) = paste(s, 'sources')
KBSC.df = reshape::melt(KBSC.df)
KBSC.df = data.frame(KBSC.df, method=rep('KeNary', length(n.control.samps)))

### Look at SVM results defined 3 ways
sim.rates = lapply(SVM.rate[3:length(FTIR.sims.results)], function(x) rowMeans(x))
# All-or-Nothing
all.or.nothing = lapply(sim.rates , function(x) ifelse(x==1, 1, 0))
all.or.nothing = lapply(all.or.nothing, mean)
all.or.nothing.df = do.call('rbind', lapply(all.or.nothing, unlist))
rownames(all.or.nothing.df) = paste(n.control.samps, 'control obs')
colnames(all.or.nothing.df) = paste(s, 'sources')
all.or.nothing.df = reshape::melt(all.or.nothing.df, id.vars=1:2)
all.or.nothing.df = data.frame(all.or.nothing.df, method=rep('SVM (all or nothing)', length(n.control.samps)))
# Per-Object
per.object = lapply(sim.rates, mean)
per.object.df = do.call('rbind', lapply(per.object, unlist))
rownames(per.object.df) = paste(n.control.samps, 'control obs')
colnames(per.object.df) = paste(s, 'sources')
per.object.df = reshape::melt(per.object.df)
per.object.df = data.frame(per.object.df, method=rep('SVM (per object)', length(n.control.samps)))
# Voting 
voting = lapply(sim.rates, function(x) ifelse(x>=2/3, 1, 0))
voting = lapply(voting, mean)
voting.df = do.call('rbind', lapply(voting, unlist))
rownames(voting.df) = paste(n.control.samps, 'control obs')
colnames(voting.df) = paste(s, 'sources')
voting.df = reshape::melt(voting.df)
voting.df = data.frame(voting.df, method=rep('SVM (voting)', length(n.control.samps)))

### Combine all results into a single data frame
all.results.df = rbind(KBSC.df, all.or.nothing.df, per.object.df, voting.df)
all.results.df = data.frame(all.results.df %>% dplyr::group_by(method, X1) %>% dplyr::mutate(Avg_Success=mean(value)))
levels(all.results.df$X2) = c(3:10)

### Create plots
ggplot(subset(all.results.df, !(X2%in%c('9 control obs', '10 control obs'))), aes(y=value, color=method)) + 
        geom_boxplot(alpha=0.5, aes(fill=method), position=position_dodge(1.)) + 
        scale_alpha_manual(values=seq(0.25, 1, length.out=8)) +
        ggthemes::scale_color_colorblind(guide="none") +
        ggthemes::scale_fill_colorblind() +
        theme(legend.position='bottom', legend.box='vertical', 
              axis.text.x=element_blank(), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
              axis.text.y=element_text(size=18), axis.title.y=element_text(size=20),
              legend.text=element_text(size=18), legend.title=element_text(size=20), title=element_text(size=21), 
              strip.text=element_text(size=18)) + 
        labs(title="Performance of KeNary vs. SVM", x="", y="Percent Correct Classification") + 
        guides(fill=guide_legend(title="Classification Method")) #+
#################################################################################################

