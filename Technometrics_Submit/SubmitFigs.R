###########################################################################################################
##### SubmitFigs.R
##### 26 November 2024
##### Create Figures in Paper
###########################################################################################################


###########################################################################################################
### Initialize Workspace 
###########################################################################################################
### Start from clean workspace 
rm(list=ls())

### Load packages 
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(dplyr)
library(data.table)
library(baseline)
library(fda)
library(rhdf5)

### Set working directory
setwd('/Users/mausdemore/Documents/GitHub/kenary/Technometrics_Submit')  

### Load data required for plots 
# Load simulated pseudo-spectra
load('Sim_Spectra.rda')
# Load data to recreate image from Ausdemore et al.
load("RMP_dat_N5_M3_fixed_trace.rData")
# Load observed paint can 19
paint_can19 = read.csv('paintcan_19.csv')
# Load observed paint can 34
paint_can34 = read.csv('paintcan_34.csv')
# Load results files 
results.files = list.files(getwd(), pattern='24.rda')
# Load source grid 
source.grid = as.matrix(data.table::fread(paste(getwd(), '/source_grid.csv', sep='')), nrow=100, ncol=10)


### Source functions
source('KeNary_Funs.R')
###########################################################################################################


###########################################################################################################
### Figure 1                                                                                              #
###########################################################################################################
# Similar spectra from two distinct paint sources. Spectra from paint can 19 in the paint dataset (black) #
#  overlaid with spectra from paint can #34 in the paint dataset (orange)                                 # 
###########################################################################################################
fig1.df = data.frame(rbind(sim.spectra[[19]], sim.spectra[[34]]))
fig1.df = cbind(PaintCan=c(rep('19', nrow(sim.spectra[[19]])), rep('34', nrow(sim.spectra[[34]]))), fig1.df)
fig1.df = reshape2::melt(fig1.df, id.vars=c(1,2))

ggplot(fig1.df, aes(x=xs.eval, y=value, col=PaintCan)) + 
 geom_line() + scale_color_colorblind() +
 ggtitle('Overlaid Spectra (Paint Can 19 versus Paint Can 34)') + 
 xlab('Wave Number') + ylab('Absorbance') +
 guides(colour=guide_legend(title="Paint Can")) + 
 theme(axis.title.x=element_text(size=17), axis.text.x=element_text(size=15), 
       axis.title.y=element_text(size=17), axis.text.y=element_text(size=15), 
       title=element_text(size=17),
       legend.position='bottom', legend.title=element_text(size=17), 
       legend.text=element_text(size=15))
###########################################################################################################



###########################################################################################################
### Figure 4a                                                                                             #
###########################################################################################################
# (a) Seven observed replicates of FTIR spectra (solid black lines) overlaid with 25 simulated replicates #
# of pseudo-spectra (solid grey lines) from paint can 19 (out of 166 paint cans comprising the dataset).  # 
# (b) Filtered spectra from paint can 19 (black) compared to filtered spectra from can 34 (orange). Areas #
# emphasized in bolded black or orange correspond to the informative areas that are considered by the     #
# kernel function given in (\ref{Eqn:KernFun}). Areas not considered by the kernel function are           #
# considered to be uninformative. By considering only the informative areas, we can better distinguish    #
# between two classes of spectra.                                                                         #
###########################################################################################################
### 4a 
spectra.19 = paint_can19
spectra.19[spectra.19==0] = min(apply(spectra.19, 2, function(x) min(x[x!=0]))[-1])
spectra.19 = cbind(xs=spectra.19[,1], -1*log10(spectra.19[,2:dim(spectra.19)[2]]/100))
spectra.19 = cbind(spectra.19[,1], 
                   t(baseline::baseline.modpolyfit(t(spectra.19[,2:ncol(spectra.19)]))$corrected))
colnames(spectra.19) = c('xs.eval', paste('obs.', 1:7, sep=''))

fig4a.orig.df = data.frame(spectra.19)
fig4a.orig.df = cbind(PaintCan=rep('Observed', nrow(spectra.19)), fig4a.orig.df)
fig4a.sim.df = data.frame(sim.spectra[[19]])
fig4a.sim.df = cbind(PaintCan=rep('Simulated', nrow(sim.spectra[[19]])), fig4a.sim.df)
fig4a.ncol = min(ncol(fig4a.orig.df)-1, ncol(fig4a.sim.df)-1)
fig4a.df = rbind(fig4a.orig.df[1:fig4a.ncol], fig4a.sim.df[1:fig4a.ncol])
fig4a.df = reshape2::melt(fig4a.df, id.vars=c(1,2))

ggplot(fig4a.df, aes(x=xs.eval, y=value, col=PaintCan)) + 
        geom_line(aes(linetype=PaintCan, linewidth=PaintCan)) + 
        scale_color_grey() + scale_linewidth_manual(values=c(1.75,0.9)) +
        ggtitle('Observed vs. Simulated Spectra (Paint Can 19)') + 
        xlab('Wave Number') + ylab('Absorbance') +
        guides(colour=guide_legend(title="Spectra Type")) + 
        theme(axis.title.x=element_text(size=17), axis.text.x=element_text(size=15), 
              axis.title.y=element_text(size=17), axis.text.y=element_text(size=15), 
              title=element_text(size=17),
              legend.position='bottom', legend.title=element_text(size=17), 
              legend.text=element_text(size=15)) + 
        labs(color="Spectra Type", linetype="Spectra Type", linewidth="Spectra Type")

### 4b
filtered.ix = index.pairs(sim.spectra[[19]][,2], sim.spectra[[34]][,2], xs.eval=sim.spectra[[19]][,1])
fig4b.df = data.frame(rbind(sim.spectra[[19]][filtered.ix,1:2], sim.spectra[[34]][filtered.ix,1:2]))
fig4b.df = cbind(PaintCan=c(rep('19', length(filtered.ix)), rep('34', length(filtered.ix))), fig4b.df)
fig4b.df = reshape2::melt(fig4b.df, id.vars=c(1,2))

ggplot(subset(fig4b.df, variable=='obs.1'), aes(x=xs.eval, y=value, col=PaintCan)) + 
        geom_line() + 
        scale_color_manual(values=c("#000000", "#E69F00"))  + 
        geom_point(data=fig4b.df, aes(x=xs.eval, y=value, col=PaintCan)) +
        ggtitle('Filtered Spectra (Paint Can 19 vs. Paint Can 34)') + 
        xlab('Wave Number') + ylab('Absorbance') +
        guides(colour=guide_legend(title="Paint Can")) + 
        theme(axis.title.x=element_text(size=17), axis.text.x=element_text(size=15), 
              axis.title.y=element_text(size=17), axis.text.y=element_text(size=15), 
              title=element_text(size=17),
              legend.position='bottom', legend.title=element_text(size=17), 
              legend.text=element_text(size=15)) + 
        labs(color="Paint Can")
###########################################################################################################



###########################################################################################################
### Figure 5                                                                                              #
###########################################################################################################
# Distribution of random match probabilities associated with a population of 166 paint cans making up the # 
# population of potential sources when we consider 5 control objects per source. Orange box plots         #
# correspond to paint cans used in the various experiments discussed in this section. See details of full #
# experiment in Ausdemore et al. 2019                                                                     #
###########################################################################################################
c.alpha5 = 0.002
source.grid = read.csv('source_grid.csv')
unique.sources = unique(unlist(source.grid))
col.ix = rep('grey40', 166)
col.ix[unique.sources] = "#E69F00"
RMP.Cluster.Results = RMP.clust
RMP = sapply(RMP.Cluster.Results, function(x,c.alpha){return(rowSums(x> c.alpha)/165)}, c.alpha=c.alpha5)
RMP = as.vector(RMP)
RMP = cbind((rep(1:166, each=20)), RMP)
colnames(RMP) = c("Spectra", "HVal")
RMP = data.frame(RMP)

ggplot(subset(RMP, Spectra<=83), aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12)) + 
        geom_boxplot(color=col.ix[1:83]) + 
        scale_color_manual(values=col.ix) +
        scale_x_discrete(breaks=pretty(rev(1:83), 83/3)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(x="", y="RMP", title="Random Match Probabilities: 5 control, 3 trace") + 
        theme(axis.title.x=element_text(size=25), axis.text.x=element_text(size=20), 
              axis.title.y=element_text(size=25), axis.text.y=element_text(size=20), 
              title=element_text(size=25),
              legend.position='bottom', legend.title=element_text(size=20), 
              legend.text=element_text(size=20)) 

ggplot(subset(RMP, Spectra>83), aes(y=HVal, x=as.factor(Spectra)))+ylim(c(0,0.12)) + 
        geom_boxplot(color=col.ix[84:166]) + 
        scale_color_manual(values=col.ix) +
        scale_x_discrete(breaks=pretty(rev(84:166), 83/3)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(x="Spectra", y="RMP", title="") + 
        theme(axis.title.x=element_text(size=25), axis.text.x=element_text(size=20), 
              axis.title.y=element_text(size=25), axis.text.y=element_text(size=20), 
              title=element_text(size=25),
              legend.position='bottom', legend.title=element_text(size=20), 
              legend.text=element_text(size=20)) 
###########################################################################################################



###########################################################################################################
### Figure 6                                                                                              #
###########################################################################################################
# Summary of KeNary performance versus SVM performance as the number of control objects increases for     #
# mi=3,4,5,6,7,8,9 control objects, and as the number of sources considered increases for n=2,3,4,5,6,7,8 # 
# sources.                                                                                                #
###########################################################################################################
### Get KBSC success rate and SVM success rate for all numbers of classes
KBSC.rate = SVM.rate = sim.rates = all.or.nothing = per.object = voting = list()
for(which.file in 1:length(results.files)){
        load(paste(getwd(), '/', results.files[which.file], sep=''))                                                        
        n.source = as.numeric(substr(paste(results.files[which.file], sep=''), 17, 17))
        
        KBSC.rate[[which.file]] = SVM.rate[[which.file]] = list()
        for(i in 3:length(FTIR.sims.results)){
                SVM.rate[[which.file]][[i]] = matrix(0, nrow=length(FTIR.sims.results[[i]]), ncol=3)
                KBSC.rate[[which.file]][[i]] = sum(unlist(lapply(FTIR.sims.results[[i]], function(x) which.max(x$BF)==x$n.sources)))/length(FTIR.sims.results[[i]])
                for(j in 1:length(FTIR.sims.results[[i]])){
                        tmp.col = apply(FTIR.sims.results[[i]][[j]]$p.SVM, 1, which.max)
                        SVM.rate[[which.file]][[i]][j,] = names(sapply(1:3, function(x, sims.results, tmp.col) 
                                sims.results[x, tmp.col[x]], FTIR.sims.results[[i]][[j]]$p.SVM, tmp.col)) == source.grid[j, n.source]
                }
        }
        
        # Look at SVM Success Rate 
        # Get per-sim success 
        sim.rates[[which.file]] = lapply(SVM.rate[[which.file]][3:length(FTIR.sims.results)], function(x) rowMeans(x))
        
        # All-or-Nothing 
        all.or.nothing[[which.file]] = lapply(sim.rates[[which.file]] , function(x) ifelse(x==1, 1, 0))
        all.or.nothing[[which.file]] = lapply(all.or.nothing[[which.file]], mean)
        
        # Per-Object 
        per.object[[which.file]] = lapply(sim.rates[[which.file]], mean)
        
        # Voting 
        voting[[which.file]] = lapply(sim.rates[[which.file]], function(x) ifelse(x>=2/3, 1, 0))
        voting[[which.file]] = lapply(voting[[which.file]], mean)
}


### Look at KBSCS for 2 classes as number of control objects varies 3:10
KBSC.df = do.call('rbind', lapply(KBSC.rate, unlist))
colnames(KBSC.df) = paste(3:10, 'control obs')
rownames(KBSC.df) = paste((1:length(results.files))+1, 'sources')
KBSC.df = reshape::melt(KBSC.df)
KBSC.df = data.frame(KBSC.df, method=rep('KeNary', length(results.files)*8))

### Look at SVM results for 2 classes as number of control objects varies 3:10
# All-or-Nothing
all.or.nothing.df = do.call('rbind', lapply(all.or.nothing, unlist))
colnames(all.or.nothing.df) = paste(3:10, 'control obs')
rownames(all.or.nothing.df) = paste((1:length(results.files))+1, 'sources')
all.or.nothing.df = reshape::melt(all.or.nothing.df, id.vars=1:2)
all.or.nothing.df = data.frame(all.or.nothing.df, method=rep('SVM (all or nothing)', length(results.files)*8))
# Per-Object
per.object.df = do.call('rbind', lapply(per.object, unlist))
colnames(per.object.df) = paste(3:10, 'control obs')
rownames(per.object.df) = paste((1:length(results.files))+1, 'sources')
per.object.df = reshape::melt(per.object.df)
per.object.df = data.frame(per.object.df, method=rep('SVM (per object)', length(results.files)*8))
# Voting 
voting.df = do.call('rbind', lapply(voting, unlist))
colnames(voting.df) = paste(3:10, 'control obs')
rownames(voting.df) = paste((1:length(results.files))+1, 'sources')
voting.df = reshape::melt(voting.df)
voting.df = data.frame(voting.df, method=rep('SVM (voting)', length(results.files)*8))

### Combine all results into a single data frame
all.results.df = rbind(KBSC.df, all.or.nothing.df, per.object.df, voting.df)
all.results.df = data.frame(all.results.df %>% dplyr::group_by(method, X1) %>% dplyr::mutate(Avg_Success=mean(value)))
levels(all.results.df$X2) = c(3:10)

### Create plots
p = ggplot(subset(all.results.df, !(X2%in%c('9 control obs', '10 control obs'))), aes(y=value, color=method)) + 
        geom_boxplot(alpha=0.5, aes(fill=method), position=position_dodge(1.)) + 
        #geom_point(aes(alpha=X2, group=method), position=position_jitterdodge(jitter.width=0.5), size=2.5) +
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

# Design block for plots 
design = c(
"
#AABBCC#
DDEEFFGG
"
)
# Create plot
p + ggh4x::facet_manual(~X1, design = design)
###########################################################################################################


###########################################################################################################
### Figure 7                                                                                              #
###########################################################################################################
# Prior and posterior probability densities for the probability of class for source combination 57, for   #
# n=2,3,5 sources, with mi=3, 8 control objects per source.                                               #
###########################################################################################################
### Define which source combination we want to consider 
which.mix = 78

### Make plots 
for(s in c(2,3,5)){
        ### Set up Experiment
        # Define number of control samples and number of test samples
        cs=3
        n.test.samps=3
        # Load data file
        load(paste(getwd(), '/', results.files[s-1], sep=''))
        # Status update
        print(paste(s, 'sources:', cs, 'control'))
        # Update number of control sampels (cs) if not a vector 
        cs = unique(cs)
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
        P.mat.control = P.mat.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], 1:sum(cs)] 
        PP.t.control = P.mat.control %*% t(P.mat.control)
        
        
        ### Prepare Adaptive MH sampler
        # Define starting point for sampler 
        theta0 = c(rep(50,choose(s,2)+s), rep(5, choose(s,2)+s), 0.1) #runif(2*(n.sources+choose(n.sources,2))+1)
        
        ### Define sources and simulation info
        which.sources = data.table::as.data.table(t(source.grid[which.mix, 1:s]))
        which.prop = paste("H", s, sep="")
        source.info = cbind(which.sources, which.prop)
        which.sources = unlist(which.sources)
        
        ### Sample spectra for control and trace objects 
        # Control Spectra
        control.spectra = NULL
        for(i in 1:s){
                control.spectra = cbind(control.spectra, 
                                         sim.spectra[[which.sources[i]]][,(1:cs[i])+1])
        }
        colnames(control.spectra) = unlist(sapply(1:s, function(x, cs){
                paste("control.", x, '.', 1:cs[x], sep="")
        }, cs))
        # Trace Spectra
        trace.spectra = sim.spectra[[which.sources[s]]][,cs[i]+1+(1:n.test.samps)]
        if(!is.matrix(trace.spectra)){
                trace.spectra = matrix(trace.spectra, ncol=1)
        }
        colnames(trace.spectra) = paste("test.", 1:n.test.samps, sep="")
        # Combine the control and the trace spectra into one matrix
        spectra.mat = cbind(control.spectra, trace.spectra)
        
        ### Get scores
        # Define spectral relevant info 
        wav.num = sim.spectra[[1]][,1]
        # At which specific points do we want to evaluate the spectra?
        xs.eval = wav.num[round(seq(1,length(wav.num),length=1000))]
        # What lag do we want to consider for the cross-correlation function?
        max.lag = 10
        # Define kernel function
        spectkern = spectral.kernel(xs.eval, max.lag, wav.num)
        # Get scores
        all.scores = apply(pairs.all, 1, function(x, spectra.mat, xs.eval, max.lag, wav.num){
                spectkern(spectra.mat[,x[1]], spectra.mat[,x[2]])
        }, spectra.mat, xs.eval, max.lag, wav.num)
        control.scores = all.scores[unlist(ix.all[1:(s+choose(s, 2))])[control.ix]]
        # Define starting point for theta0
        n.samps=100000
        theta0.mat=matrix(rep(theta0, n.samps), nrow=n.samps, 
                          ncol=2*(s+choose(s,2))+1, byrow=T) + 
                matrix(c(runif(2*(s+choose(s,2))*n.samps,-1,1), runif(1*n.samps,0,1)), 
                       ncol=2*(s+choose(s,2))+1, nrow=n.samps)
        
        ### Create Prior Probability Plots
        BF.eval.time = system.time({
                BF.samps = sapply(1:s, BF.eval, theta0.mat, all.scores, P.mat.all, s, cs, n.test.samps)
        })
        tmp.df = exp(BF.samps/100)
        tmp.rowsum = rowSums(tmp.df)
        for(i in 1:s){
                tmp.df[,i] = tmp.df[,i]/tmp.rowsum
        }
        tmp = reshape2::melt(tmp.df)[,-1]
        tmp.df = data.frame(tmp, control=rep(ncol(tmp.df), nrow(tmp)))
        
        p1 = ggplot(subset(tmp.df), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
                stat_density(alpha=0.5, trim=T) + 
                scale_fill_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")[1:s]) + 
                scale_color_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")[1:s], guide="none") + 
                guides(fill=guide_legend(title="Probability of Class... "), col=guide_legend(title="Probability of Class... ")) + 
                labs(title=paste("Prior Probability for Combination ", which.mix, ' (', s, ' sources)', sep='')) + 
                xlab("") +
                theme(legend.position='bottom', legend.box='vertical', 
                      axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
                      axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
                      legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
                      strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
                facet_wrap(~control, scales="free", nrow=1) + 
                scale_x_continuous(n.breaks=5) 
        
        ### Create Posterior Probability Plots
        tmp=list()
        for(i in 3:10){
                tmp.df = exp(FTIR.sims.results[[i]][[which.mix]]$BF.samps/c(10000, 25000, 150000)[ifelse(s==2 | s==3, s-1, 3)])
                tmp.rowsum = rowSums(tmp.df)
                
                for(j in 1:ncol(tmp.df)){
                        tmp.df[,j] = tmp.df[,j]/tmp.rowsum
                }
                tmp[[i]] = reshape2::melt(tmp.df)[,-1]
                tmp[[i]] = cbind(tmp[[i]], control=rep(i, nrow(tmp[[i]])))
        }
        tmp.df = do.call('rbind', tmp)
        
        p2 = ggplot(subset(tmp.df, control%in%c(3,8)), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
                stat_density(alpha=0.5) + 
                scale_fill_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")[1:s])  + 
                scale_color_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")[1:s] , guide="none") + 
                guides(fill=guide_legend(title="Probability of Class... "), col=guide_legend(title="Probability of Class... ")) + 
                labs(title=paste("Distribution of Posterior Probability of Class for Combination ", which.mix, ' (', s, ' sources)', sep='')) + 
                xlab("") +
                theme(legend.position='bottom', legend.box='vertical', 
                      axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
                      axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
                      legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
                      strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
                facet_wrap(~control, scales="free", nrow=1) + 
                scale_x_continuous(n.breaks=3) 
        
        grid.arrange(p1, p2, nrow=1, widths=c(0.75, 1.5))
}
###########################################################################################################


###########################################################################################################
### Figure 9                                                                                              #
###########################################################################################################
# SP-AMS spectra of soot recovered from detonation in air (black), from detonation in Argon (yellow), and #
# from a Fullerene standard (blue).                                                                       #
###########################################################################################################

### Load hdf5 files
# Air Spectra
RawSpectra_Air_tmp <- h5read(file=paste(getwd(), '/CompB_Air_1d_V_BC.h5', sep=''), name='Raw_Spectra_Diff/TodoRunRawSpectra')
# Ar Spectar
RawSpectra_Ar_tmp <- h5read(file=paste(getwd(), '/CompB_Ar_1d_V_BC.h5', sep=''), name='Raw_Spectra_Diff/TodoRunRawSpectra')
# Fullerene Spectra
RawSpectra_Ful_tmp <- h5read(file=paste(getwd(), '/CompB_Ful_1d_V_BC.h5', sep=''), name='Raw_Spectra/TodoRunRawSpectraDiff')
# Unknown Spectra
Unknown_Data_tmp <- h5read(file=paste(getwd(), '/D160707.h5', sep=''), name='Raw_Spectra_Diff/TodoRunRawSpectra')

### Get union of time stamps
time_union <- union(union(union(RawSpectra_Ar_tmp[1,1,], RawSpectra_Air_tmp[1,1,]), RawSpectra_Ful_tmp[1,1,]), Unknown_Data_tmp[1,1,])
time_union <- time_union[order(time_union)]
ArTime_ix <- which(time_union %in% RawSpectra_Ar_tmp[1,1,])
AirTime_ix <- which(time_union %in% RawSpectra_Air_tmp[1,1,])
FulTime_ix <- which(time_union %in% RawSpectra_Ful_tmp[1,1,])
UnknownTime_ix <- which(time_union %in% Unknown_Data_tmp[1,1,])
RawSpectra_Ar <- matrix(0, nrow=length(time_union), ncol=dim(RawSpectra_Ar_tmp)[1]+1)
RawSpectra_Air <- matrix(0, nrow=length(time_union), ncol=dim(RawSpectra_Air_tmp)[1]+1)
RawSpectra_Ful <- matrix(0, nrow=length(time_union), ncol=dim(RawSpectra_Ful_tmp)[1]+1)
RawSpectra_Unk <- matrix(0, nrow=length(time_union), ncol=dim(Unknown_Data_tmp)[1]+1)
RawSpectra_Ar[,1] <- time_union
RawSpectra_Air[,1] <- time_union
RawSpectra_Ful[,1] <- time_union
RawSpectra_Unk[,1] <- time_union

### Organize spectra into a matrix
RawSpectra_Air[AirTime_ix, 2:(dim(RawSpectra_Air_tmp)[1]+1)] <- t(RawSpectra_Air_tmp[1:dim(RawSpectra_Air_tmp)[1],2,])
RawSpectra_Ar[ArTime_ix, 2:(dim(RawSpectra_Ar_tmp)[1]+1)] <- t(RawSpectra_Ar_tmp[1:dim(RawSpectra_Ar_tmp)[1],2,])
RawSpectra_Ful[FulTime_ix, 2:(dim(RawSpectra_Ful_tmp)[1]+1)] <- t(RawSpectra_Ful_tmp[1:dim(RawSpectra_Ful_tmp)[1],2,])
RawSpectra_Unk[UnknownTime_ix, 2:(dim(Unknown_Data_tmp)[1]+1)] <- t(Unknown_Data_tmp[1:dim(Unknown_Data_tmp)[1],2,])

### Create plot
tmp.df = rbind(reshape::melt(data.frame(RawSpectra_Air), id.vars=1), 
               reshape::melt(data.frame(RawSpectra_Ar), id.vars=1),
               reshape::melt(data.frame(RawSpectra_Ful), id.vars=1),
               reshape::melt(data.frame(RawSpectra_Unk), id.vars=1))
tmp.df = cbind(tmp.df, Sample=rep(c("Air", "Argon", "Fulerene", "Unknown"), times=c(2575000, 2832500, 4281680, 2317500)))

ggplot(subset(tmp.df, Sample!="Unknown" & variable=="X2"), aes(x=X1, y=log(value+50), col=Sample)) + 
        geom_line(linewidth=0.5) +
        ylab("log(Intensity + 50)") + xlab("m/z") + ggtitle("SP-AMS Spectra") + 
        ggforce::facet_zoom(xlim=c(0,250)) + 
        ggthemes::scale_color_colorblind() +
        theme(legend.position="bottom", 
              legend.text=element_text(size=15), 
              legend.title=element_text(size=20), 
              axis.text.x=element_text(size=15), 
              axis.text.y=element_text(size=15), 
              strip.text=element_text(size=17), 
              title=element_text(size=20)) 
###########################################################################################################


###########################################################################################################
### Figure 10                                                                                             #
###########################################################################################################
# Posterior probability of detonation class associated with the unknown set of trace samples.             # 
###########################################################################################################

### Define spectral relevant info 
# How many basis functions do we want to use to define the spectra?
p.basis <- 300
# How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
# At which specific points do we want to evaluate the spectra?
wav.num <- time_union
xs.eval <- wav.num[round(seq(1,length(wav.num),length=p.eval))]
# What lag do we want to consider for the cross-correlation function?
max.lag <- 10
sp.ams.dat <- list(RawSpectra_Air[,-1], RawSpectra_Ar[,-1], RawSpectra_Ful[,-1], RawSpectra_Unk[,-1])

### Status update
s <- 3
cs <- c(ncol(RawSpectra_Air), ncol(RawSpectra_Ar), ncol(RawSpectra_Ful))-1
n.test.samps = ncol(RawSpectra_Unk)-1
### Update cs.samps if not a vector 
cs <-unique(cs)
if(length(cs)==1){
        cs <- rep(cs, s)
}

### Define object pairwise comparisons for parameter estimation
# Get an initial grid of pairwise comparisons 
pairs.all <- t(combn(1:(sum(cs)+n.test.samps),2))
# Define indices to distinguish between trace and control objects - which.prop here doesn't matter 
ix.all <- get.indices(n.sources=s, n.sizes=cs, 
                      n.unknown=n.test.samps, which.prop=paste("H", s, sep=''))
# Grab control pairs 
control.ix <- which(unlist(ix.all[1:(s+choose(s, 2))])%in% ix.all$control.IX)
control.pairs <- pairs.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], ]


### Create P matrices
P.mat.all <- make.P(sum(cs)+n.test.samps)
P.mat.control <- P.mat.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], 1:sum(cs)] 
PP.t.control <- P.mat.control %*% t(P.mat.control)



### Prepare Adaptive MH sampler
# Define starting point for sampler 
theta0 <- c(rep(50,choose(s,2)+s), rep(5, choose(s,2)+s), 0.1) #runif(2*(s+choose(s,2))+1)
kernfun <- Metrics::mse
sims.info=c(1,2,3)
n.samps=5000*5
seed=NULL
SVM.analysis=FALSE

### Define sources and simulation info
which.sources <- (t(sims.info[1:s]))
source.info <- cbind(which.sources)

### Sample spectra for control and trace objects 
# sample cs.samps spectra for each control source
if(!is.null(seed)) set.seed(seed)
control.spectra <- NULL
for(i in 1:s){
        control.spectra <- cbind(control.spectra, 
                                 sp.ams.dat[[which.sources[i]]])
}
colnames(control.spectra) <- unlist(sapply(1:s, function(x, cs){
        paste("control.", x, '.', 1:cs[x], sep="")
}, cs))
# Sample n.test.samps spectra for the trace source
trace.spectra <- sp.ams.dat[[s+1]]
if(!is.matrix(trace.spectra)){
        trace.spectra <- matrix(trace.spectra, ncol=n.test)
}
colnames(trace.spectra) <- paste("test.", 1:n.test.samps, sep="")
# Combine the control and the trace spectra into one matrix
spectra.mat <- cbind(control.spectra, trace.spectra)


### Get scores
all.scores <- apply(pairs.all, 1, function(x, spectra.mat){
        kernfun(spectra.mat[,x[1]], spectra.mat[,x[2]])
}, spectra.mat)
all.scores=log(all.scores)
control.scores <- all.scores[unlist(ix.all[1:(s+choose(s, 2))])[control.ix]]


### Get Model Parameters via Adaptive MH
# Define function for calculating log-posterior
MCMC.logpost <- logpost.scores(control.scores, n.sources=s, n.samps=cs, PP.t=PP.t.control) 
param.samp.time <- system.time({
        param.samps <- MCMC(MCMC.logpost, init=theta0, adapt=TRUE, 
                            n=n.samps, acc.rate=0.234, 
                            scale=rep(1, 2*(s+choose(s,2))+1))$samples
})
param.samps.tail <- param.samps[seq(0.75*nrow(param.samps), nrow(param.samps), 5),]

### Evaluate Bayes Factor Using Sampled Parameter Values 
BF.eval.time <- system.time({
        BF.samps <- sapply(1:s, BF.eval, param.samps.tail, all.scores, P.mat.all, s, 
                           cs, n.test.samps)
})
BF.num <- apply(BF.samps, 2, function(x) mean(exp(x[is.finite(x)]/100))) 
BF <- sapply(BF.num, function(x, BF.num) x/sum(BF.num), BF.num)

### Get ready to plot
tmp.df = exp(BF.samps/1000)
tmp.rowsum = rowSums(tmp.df)
for(i in 1:3){
        tmp.df[,i] = tmp.df[,i]/tmp.rowsum    
}
tmp = reshape2::melt(tmp.df)[,-1]
tmp.df = as.data.frame(tmp)

### Create plot
ggplot(tmp.df, aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
        stat_density(alpha=0.5) + 
        ggthemes::scale_fill_colorblind(labels=c("Detonation in Air", "Detonation in Argon", "Fullerene Standard")) + 
        ggthemes::scale_color_colorblind(guide="none") + 
        guides(fill=guide_legend(title="Probability of... ")) + 
        labs(title=paste("Distribution of Probability of Class for Unknown Soot Sample")) + 
        xlab("") +
        theme(legend.position='bottom', legend.box='vertical', 
              axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
              axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
              legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
              strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
        scale_x_continuous(n.breaks=6)
###########################################################################################################
