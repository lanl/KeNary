#################################################################################################
##### PaintSims_Study.R
##### 26 November 2024
##### Simulations to Run on a Cluster
#################################################################################################


#################################################################################################
##### Initialize Workspace 
#################################################################################################
### Start from a clean workspace 
rm(list=ls())

#################################################################################################
###                                       Info to Change                                        #
#################################################################################################
### Define paths                                                                                #
# Set working directory - this is the folder where your data files are saved                    #
setwd('/Users/mausdemore/Documents/GitHub/kenary/Technometrics_Submit')                         #                                                                      #
# Source your function file                                                                     #
source("KeNary_Funs.R")                                                                         #
                                                                                                #
### Define the number of cores that should be used                                              #
n.cores <- parallel::detectCores()-2                                                            #
                                                                                                #
### Define the parameters for the study                                                         #
# Define the number of sources to be studied                                                    #
n.sources <-  7                                                                                 #
# Define the number of control objects - this can be a vector to study multiple sample sizes    #
n.control.samps <- 10                                                                           #
# Define the number of test objects - this is the set of objects whose source is unknown        #
n.test.samps <- 3                                                                               #
#################################################################################################
###                       !!! Do not change anything below this line !!!                        #
#################################################################################################



#################################################################################################
##### Load the Data 
#################################################################################################
### Load paint data 
load('Sim_Spectra.rda')
# Get wave numbers from data (x-axis)
wav.num <- sim.spectra[[1]][,1]

### Load source information
source.grid <- as.matrix(fread('source_grid.csv'), nrow=100, ncol=10)
#################################################################################################

#################################################################################################
##### Set up Model 
#################################################################################################
### Define relevant spectral info
# How many basis functions do we want to use to define the spectra?
p.basis <- 300
# How many points do we want to use to evaluate the pseudo-spectra?
p.eval <- 1000
# At which specific points do we want to evaluate the spectra?
xs.eval <- wav.num[round(seq(1,length(wav.num),length=p.eval))]
# What lag do we want to consider for the cross-correlation function?
max.lag <- 10

### Define kernel function for spectra
spectkern <- spectral.kernel(xs.eval, max.lag, wav.num)

### Define list to store data 
FTIR.sims.results <- list()

for(cs in n.control.samps){
     ### Status update
     print(cs)
     s <- n.sources
     
     ### Update n.control.samps if not a vector 
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
     P.mat.control <- P.mat.all[unlist(ix.all[1:(s+choose(s, 2))])[control.ix], 
                                1:sum(cs)] 
     PP.t.control <- P.mat.control %*% t(P.mat.control)
     
     
     ### Prepare Adaptive MH sampler
     # Define starting point for sampler 
     theta0 <- c(rep(50,choose(s,2)+s), rep(5, choose(s,2)+s), 0.1) 
     
     ### Get performance of Adaptive MH versus SVM for paint data
     registerDoParallel(cores=n.cores)
     FTIR.sims.results <- foreach(i=1:nrow(source.grid)) %dopar%
                                FTIR.sims(sims.info=source.grid[i,], dat.ls=sim.spectra, n.sources=s, n.control=cs,
                                          n.test=n.test.samps, control.pairs=control.pairs, 
                                          PP.t.control=PP.t.control, control.ix=control.ix, P.mat.all=P.mat.all, 
                                          theta0=theta0, package.vect=package.vect, kernfun=spectkern, wav.num=wav.num, 
                                          SVM.analysis=TRUE, seed=0806, n.samps=5e05)
     stopImplicitCluster()
     save(FTIR.sims.results, file=paste('FTIRSimsResults_',
                                        n.sources, 'sources_',
                                        unique(cs), 'control_',
                                        n.test.samps, 'test_',
                                        stringr::str_replace_all(format(Sys.time(), "%e%b%y"), " ", ''),
                                        '.rda', sep=''))
}
#################################################################################################

