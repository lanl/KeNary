#################################################################################################
##### KeNary_Funs.R
##### 26 November 2024
##### Define Functions for Simulations
#################################################################################################


#################################################################################################
##### Load necessary libraries 
#################################################################################################
### Define packages
package.vect <- c('adaptMCMC', 'LaplacesDemon', 'kernlab', 'baseline', 'fda', 'Matrix', 
                  'mvtnorm', 'mvnfast', 'here', 'data.table', 'bayesmeta', 'doParallel')

# If package is not installed: install; Load 
package.check <- lapply(package.vect, function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x, dependencies = TRUE)
  }
 library(x, character.only = TRUE)
})
#################################################################################################

 
#################################################################################################
### Function:  logpost.scores
#################################################################################################
### * Output:  log-posterior of score vector 
### * Inputs: 
###            y:         vector of scores: 
###            n.sources: number of sources being considered; 
###            n.samps:   number of samples per source; 
###            PP.t:      design matrix for constructing covariance matrix
#################################################################################################
logpost.scores = function(y, n.sources, n.samps, PP.t){
  x.fct = function(par){
    ### Adjust theta parameters to be greater than zero if sampled less than zero
    theta = par[1:(n.sources+choose(n.sources,2))]
    theta <- abs(theta)
    
    ### Adjust sigma params to be greater than zero if sampled less than zero
    sigs = par[(choose(n.sources,2)+n.sources)+1:(n.sources+choose(n.sources,2))]
    sigs <- abs(sigs)

    ### Adjust sig.e term to be between zero and one if sampled outside zero and one
    sig.e = par[length(par)]
    sig.e = invlogit(sig.e)
    sig.e = ifelse(sig.e==1, 1-1e-5, sig.e)
    sig.e = ifelse(sig.e==0, 1e-5, sig.e)
    
    ### Do matrix algebra for Delta, Sigma and Theta
    Delta <- diag(rep(sigs, times=c(choose(n.samps,2), apply(combn(n.samps,2),2,prod))))
    Sigma.tmp = sweep((PP.t*(1-sig.e)/2 + diag(sig.e, nrow(PP.t))),diag(Delta),MARGIN=2,FUN="*")
    Sigma = sweep(Sigma.tmp,diag(Delta),MARGIN=1,FUN="*")
    Theta = rep(theta, times=c(choose(n.samps,2), apply(combn(n.samps,2),2,prod)))
    # Get Loglikelihood
    mvnfast::dmvn(y, mu=Theta,sigma=base::chol(Sigma), log=TRUE, isChol=T) + 
      # log prior on sig.e
      dbeta(sig.e, 0.5, 0.5, log=TRUE) + 
      # lomax prior on theta
      sum(dlomax(theta, shape=6, scale=1)) +
      # halfcauchy prior on sigs
      sum(dhalfcauchy(sigs, 1e-5, log=TRUE))
  }
  return(x.fct)
}
#################################################################################################


#################################################################################################
### Function:   get.indices
#################################################################################################
### * Output:  indices associated with different scoring groups 
### * Inputs:
###            n.sources:  number of sources considered; 
###            n.sizes:    number of objects per source;
###            n.unknown:  number of test objects;
###            which.prop: utility of function ("param.est" for parameter estimation; "which.prop"  
###                        for BF evaluation);
#################################################################################################
get.indices <- function(n.sources, n.sizes, n.unknown=0, which.prop=c("param.est", which.prop)){
  ### NOTE:
  ### The trace objects are ALWAYS the last objects to be compared!
  
  ### Define a grid of pair-wise object comparisons
  my.grid <- t(combn(n.unknown+sum(n.sizes),2))
  
  ### If we are doing parameter estimation, consider the scenario in which we have no
  ### trace objects
  if (which.prop == "param.est"){
    ### When n.eu == 0, proceed
    if (n.unknown==0){
      ### Create a simple vector of the endpoints of the control objects
      n.samples.vect <- cumsum(c(n.sizes))
      
      ### Get the indices for objects that correspond to each considered source
      source.ls <- lapply(1:length(n.samples.vect), function(x, n.samples.vect, n.sizes){
        (n.samples.vect[x]-(n.sizes[x]-1)):n.samples.vect[x]}, n.samples.vect, n.sizes)
      
      ### Get the indices that correpond to each within-source comparison scenarios
      WI.IX <- lapply(source.ls, function(x, my.grid){
        which(my.grid[,1] %in% x & my.grid[,2] %in% x)
      }, my.grid)
      
      ### Get the indices that correspond to each of the between-source comparison scenarios
      BT.mat.ix <- combn(n.sources,2)
      BT.IX <- unlist(apply(BT.mat.ix, 2, function(x, my.grid, source.ls){
        list(which((my.grid[,1] %in% source.ls[[x[1]]] & my.grid[,2] %in% source.ls[[x[2]]]) | (my.grid[,1] %in% source.ls[[x[2]]] & my.grid[,2] %in% source.ls[[x[1]]])))
      }, my.grid, source.ls), recursive=FALSE)
      
      ### Identify which objects are trace objects and which are control objects
      trace.IX <- NULL
      control.IX <- (1:choose(sum(n.sizes),2))
    }
    ### If n.unknown != 0, warn!
    else {
      stop("STOP! Number of trace objects needs to be 0")
    }
  } 
  
  ### If we are considering indices under a certain hypothesis, consider the scenario in which
  ### we have n.eu trace objects in addition to n0 control objects per source
  else{
    ### Define the source that the trace objects originate from
    n.tmp <- as.numeric(substr(which.prop, 2, nchar(which.prop)))
    
    ### Create a simple vector of the endpoints of the control objects
    n.samples.vect <- cumsum(c(n.sizes))
    
    ### Get the indices for objects that correspond to each considered source (including n.eu)
    source.ls <- lapply(1:length(n.samples.vect), function(x, n.samples.vect, n.sizes){
      (n.samples.vect[x]-(n.sizes[x]-1)):n.samples.vect[x]
    }, n.samples.vect, n.sizes)
    source.ls[[n.tmp]] <- c(source.ls[[n.tmp]], 
                            (n.samples.vect[n.sources]+1):(n.samples.vect[n.sources]+n.unknown))
    
    ### Get the indices that correpond to each within-source comparison scenarios
    WI.IX <- lapply(source.ls, function(x, my.grid){
      which(my.grid[,1] %in% x & my.grid[,2] %in% x)
    }, my.grid)
    
    ### Get the indices that correspond to each of the between-source comparison scenarios
    BT.mat.ix <- combn(n.sources,2)
    BT.IX <- apply(BT.mat.ix, 2, function(x, my.grid, source.ls){
      which((my.grid[,1] %in% source.ls[[x[1]]] & my.grid[,2] %in% source.ls[[x[2]]]) | 
              (my.grid[,1] %in% source.ls[[x[2]]] & my.grid[,2] %in% source.ls[[x[1]]]))
    }, my.grid, source.ls)
    if(ncol(BT.mat.ix)==1){
      BT.IX <- list(c(BT.IX))
    }
    
    ### Identify which objects are trace objects and which are control objects
    trace.IX <- which(my.grid[,1] > n.samples.vect[n.sources] | my.grid[,2] > n.samples.vect[n.sources])
    control.IX <- (1:dim(my.grid)[1])[-trace.IX]
  }
  
  
  ### Return list of various indices
  return(c(WI.IX=WI.IX, BT.IX=BT.IX, list(trace.IX=trace.IX, control.IX=control.IX)))
}
#################################################################################################


#################################################################################################
### Function:  sample.FTIR.objects  
#################################################################################################
### * Output:  additional spectra sampled using B-spline basis functions 
### * Inputs: 
###            p.basis: the number of basis functions used to create the pseudo-spectra; 
###            p.eval:  the number of points at which the spectra will be evaluated 
###                     (default = 500);
###            spectra: the matrix of training spectra from a single source where the first   
###                     column contains values for the x-axis; 
###            n.obs:   the number of objects to sample from the source of the matrix  
###                     of spectra;
###            p.pca:   number of pca used to define the spectra (default = 7);
###            xs.eval: the vector of points at which the pseudo-spectra will be evaluated 
###                     (default = NULL);
#################################################################################################                                
sample.FTIR.objects <- function(p.basis, p.eval=500, spectra, n.obs, p.pca=7, xs.eval=NULL){
  ### Get the range of the x-axis 
  range.spectra <- range(spectra[,1])
  
  ### Sort the spectra by ascending x-axis
  ix.spectra <- order(spectra[,1])
  spectra <- spectra[ix.spectra,]
  
  ### Get evaluation points by ascending order
  if (!is.null(xs.eval)){
    ix.spectra <- ix.spectra[1:length(xs.eval)]-min(ix.spectra[1:length(xs.eval)])+1
    xs.eval <- xs.eval[ix.spectra]
  } else {
    ix.spectra <- ix.spectra[1:p.eval]-min(ix.spectra[1:p.eval])+1
    xs.eval <- seq(spectra[1,1],spectra[dim(spectra)[1],1],length=p.eval)
  }
  
  ### Model basis functions using 'fda' package
  # create the basis functions
  bspline.basis <- create.bspline.basis(rangeval=range(spectra[,1]),nbasis=p.basis, norder=4)
  # create the basis evaluation points
  basis.eval.points <- seq(range(spectra[,1])[1], range(spectra[,1])[2],length=dim(spectra)[1])
  # evaluate the basis functions at evaluation points
  basis.evals <- eval.basis(basis.eval.points, bspline.basis)
  # get the basis functions coefficients 
  basis.coef <- solve(crossprod(basis.evals), 
                      crossprod(basis.evals, data.matrix(spectra[,2:dim(spectra)[2]])))
  # create functional objects
  spectra.fd <- fd(basis.coef, bspline.basis)
  
  ### PCA Weighting 
  # functional pca decomp 
  pca.spectra <- pca.fd(spectra.fd, p.pca)
  # get weighted coefficients 
  weights <- rep(0, p.pca)
  # reconstruct coefficients using pca weights
  reconstructed.coef <- rowMeans((pca.spectra$harmonics$coefs)%*%(t(pca.spectra$scores)+weights)) + 
    pca.spectra$meanfd$coefs
  
  ### Sample spectra using basis functions from above
  # define mean and covariance for sampling based on observed basis.coefs
  mean.par <- rowMeans(basis.coef)
  cov.par <- (cov(t(basis.coef)))
  if (sum(diag(cov.par))<0.001){
    cov.par <- cov.par*0.001/sum(diag(cov.par))
  }
  
  ### Sample coefficients for objects 
  basiscoef.sample <- t(rmvnorm(n.obs, mean.par, cov.par))
  
  ### Create the functional objects
  tmp.fd <- fd(basiscoef.sample, bspline.basis)
  tmp.spectra.mat <- cbind(xs.eval, as.matrix((eval.fd(xs.eval, tmp.fd))))
  tmp.spectra.mat <- tmp.spectra.mat[ix.spectra,]
  
  ### Name objects
  obs.names <- paste('obs', 1:n.obs, sep=".")
  colnames(tmp.spectra.mat) <- c("xs.eval", obs.names)
  
  ### Return sampled objects 
  return(tmp.spectra.mat)
}
#################################################################################################                                


#################################################################################################                                
### Function:  make.P
#################################################################################################
### * Output:  P matrix of dimension choose(n.obs, 2)
### * Inputs: 
###            n.obs: number of objects to be compared; 
#################################################################################################                                
make.P <- function(n.obs){
  n.exp <- t(combn(n.obs,2))
  P <- apply(n.exp, 1, function(x,n.obs){tmp <- rep(0,n.obs); tmp[x] <- 1; return(tmp)}, n.obs)
  return(t(P))
}
#################################################################################################                                


#################################################################################################                                
### Function:  nonstat.kern
################################################################################################# 
### * Output:  comparison to two objects via a non-stationary kernel (square-root of product)
### * Inputs:  
###            X1: first object to be compared;
###            X2: second object to be compared;
################################################################################################# 
nonstat.kern <- function(X1, X2){
  return(sqrt(X1*X2))
}
#################################################################################################


#################################################################################################                                
### Function:  poly.kern
################################################################################################# 
### * Output:  comparison to two objects via a non-stationary kernel (square-root of product)
### * Inputs:  
###            X1: first object to be compared;
###            X2: second object to be compared;
################################################################################################# 
poly.kern <- polydot(degree=1)
#################################################################################################


################################################################################################# 
### Function:  spectral.kernel
#################################################################################################
### * Output:  comparison of two FTIR spectra via a self-defined kernel
### * Inputs: 
###            xs.eval: xs at which spectra are evaluated;
###            max.lag: maximum lag for cross-correlation function;
###            wav.num: wave numbers at which spectra are to be evaluated;
################################################################################################# 
spectral.kernel <- function(xs.eval, max.lag, wav.num){
     spec.val <- function(X1, X2){
          # Get overlapping indices 
          ix.pair <- index.pairs(X1, X2, xs.eval=xs.eval)
          
          # Get cross-correlation
          tmp.ccf <- ccf(X1[ix.pair], X2[ix.pair], lag.max=max.lag, plot=FALSE)$acf
          tmp.ccf <- (2*max.lag+1 - sqrt(tmp.ccf%*%tmp.ccf))*10
          
          # Get Euclidean Norm
          X1.norm <- (X1[ix.pair]-min(X1[ix.pair]))
          X1.norm <- X1.norm/(sum(X1.norm*mean(abs(diff(wav.num)))))*length(X1.norm)
          X2.norm <- (X2[ix.pair]-min(X2[ix.pair]))
          X2.norm <- X2.norm/(sum(X2.norm*mean(abs(diff(wav.num)))))*length(X2.norm)
          tmp.euclid.norm <- sqrt(t(X1.norm-X2.norm)%*%(X1.norm-X2.norm))/length(X1.norm)*10000
          
          spec.score <- log(tmp.ccf*tmp.euclid.norm)*10
          
          return(ifelse(is.finite(spec.score), spec.score, 0))    
     }
     return(new('kernel', .Data=spec.val, kpar=list(xs.eval=xs.eval, max.lag=max.lag, wav.num=wav.num)))
}
################################################################################################# 


################################################################################################# 
### Function:  index.pairs
#################################################################################################
### * Output:  indices of a pair of spectra to be considered 
### * Inputs: 
###            X1:      first baseline corrected spectra to be compared; 
###            X2:      second baseline corrected spectra to be compared;
###            xs.eval: xs at which spectra are evaluated;
################################################################################################# 
index.pairs <- function(X1, X2, xs.eval){
  # Get indices to keep for first spectra
  ix.X1 <- index.indivs(X1, xs.eval)
  # Get indices to keep for second spectra
  ix.X2 <- index.indivs(X2, xs.eval)
  # Get union of two sets of indices
  ix <- sort(union(ix.X1, ix.X2))
  
  return(ix)
}
################################################################################################# 


################################################################################################# 
### Function:  index.indiv
################################################################################################# 
### * Output:  indices of a single spectra to consider in comparison
### * Inputs: 
###            X:       a baseline corrected spectra;
###            wav.num: wave numbers at which spectrum is to be evaluated;
################################################################################################# 
index.indivs <- function(X, wav.num){
  # Get information about X 
  n <- round(length(X)/50)
  quantile.param <- 0.6
  
  # Get first derivate
  first.deriv <- diff(X)/mean(diff(wav.num))*100
  
  # Get moving average for intensity (we keep it large to 'protect' the peaks)
  moving.avg.intensity <- as.vector(stats::filter(X, rep(1/n,n)))
  moving.avg.intensity[is.na(moving.avg.intensity)] <- 0
  
  # Get moving average for derivative (we keep it large to 'protect' the peaks)
  moving.avg.deriv <- as.vector(stats::filter(c(0, first.deriv),rep(1/n, n)))
  moving.avg.deriv[is.na(moving.avg.deriv)] <- 0
  
  # Get epsilon for intensity threshold
  epsilon.intensity <- quantile(moving.avg.intensity, quantile.param, na.rm=TRUE)
  # Get epsilon for derivative threshold
  epsilon.deriv <- quantile(abs(moving.avg.deriv), quantile.param, na.rm=TRUE)
  
  # Get indices of points we want to keep for intensity
  ix.intensity <- which(moving.avg.intensity>epsilon.intensity)
  # Get indices of points we want to keep for derivative
  ix.deriv <- which(abs(moving.avg.deriv)>epsilon.deriv)
  
  # Get indices we keep from both
  ix.keep <- sort(union(ix.deriv, ix.intensity))
  
  return(ix.keep)
}
#################################################################################################

################################################################################################# 
### Function:  BF.eval
################################################################################################# 
### * Output:  Bayes Factor evaluation 
### * Inputs:
###            which.prop:      proposition under which Bayes Factor is to be evaluated;
###            param.samps:     posterior samples from adaptive MH step;
###            all.scores:      vector of all pairwise comparisons between all objects;
###            P.mat.all:       original P matrix for all object comparisons;
###            n.sources:       number of sources being considered 
###            n.control.samps: number of control samples considered, per sources;
###            n.test.samps:    number of test samples considered;
################################################################################################# 
BF.eval <- function(prop.num, param.samps, all.scores, P.mat.all,
                    n.sources, n.control.samps, n.test.samps){
  
  # define indices to put grid in correct order 
  test.ix <- get.indices(n.sources=n.sources, 
                         n.sizes=n.control.samps, 
                         n.unknown=n.test.samps, 
                         which.prop=paste('H', prop.num, sep=''))
  
  # get scores (trace scores at the end)
  test.scores <- all.scores[unlist(test.ix[1:(n.sources+choose(n.sources,2))])]
  
  # make P matrix 
  P.mat.test <- P.mat.all[unlist(test.ix[1:(n.sources+choose(n.sources,2))]),]
  PP.t.test <- P.mat.test %*% t(P.mat.test)
  
  # get log-posterior for all samples 
  all.samps <- n.control.samps
  all.samps[prop.num] <- n.control.samps[prop.num] + n.test.samps
  test.logpost <- logpost.scores(test.scores, 
                                 n.sources=n.sources, 
                                 n.samps=all.samps, 
                                 PP.t.test)
  HX <- apply(param.samps, 1, test.logpost)
}
################################################################################################# 


################################################################################################# 
### Function: FTIR.sims 
################################################################################################# 
### * Output: list of sampled parameters (param.samps), samples used in BF calculation (BF.samps), 
###           initial samples for inital Bayes Factor (BF.samps.init), numerator of Bates Factor 
###           (BF.num), calculated Bayes Factor (BF), numerator of initial Bayes Factor 
###           (BF.num.init), calculated inital Bayes Factor (BF.init), number of sources 
###           considered (n.sources), time to sample parameters (param.samp.time), time to  
###           evaluate Bayes Factor (BF.eval.time)
### * Input: 
###           sims.info:     vector of information for simulations, consisting of sources to be 
###                          compared, and which Hypothesis is "true"
###           dat.ls:        list of spectra to pull from; each list object is a different   
###                          source of paint in the form of a matrix where rows correspond to    
###                          wave numbers and columns correspond to different samples of paint
###           n.sources:     number of sources to be considered
###           n.control:     number of control objects per source
###           n.test:        number of test objects observed
###           control.pairs: matrix of pairs of objects to be compared
###           PP.t.control:  PP.t matrix associated with control objects 
###           control.ix:    indices corresponding to control objects 
###           P.mat.all:     P matrix for all objects
###           theta0:        initial parameter vector for Adaptive MH algorithm
###           package.vect:  packages used in parallel computing
###           kernfun:       kernel function to be used to compare objects
###           n.samps:       number of samples to be obtained from Adaptive MH sampler 
###                          (default=100000)
###           wav.num:       wave numbers to be considered
###           max.lag:       maximum lag to consider in spectkern (default=10)
###           p.basis:       the number of basis functions used to create the pseudo-spectra 
###                          (default = 300); 
###           p.eval:        the number of points at which the spectra will be evaluated 
###                          (default = 1000);
###           SVM.analysis:  should SVM analysis be perfromed? (default=FALSE)
###           seed:          do we want to set a seed? If so, seed=XXX where XXX is the seed to 
###                          be used to replicate experiments (default=NULL)
################################################################################################# 
FTIR.sims <- function(sims.info, dat.ls, n.sources, n.control, n.test,
                      control.pairs, PP.t.control, control.ix, P.mat.all,
                      theta0, package.vect, kernfun, n.samps=100000,
                      wav.num, max.lag=10, p.basis=300, p.eval=1000, 
                      SVM.analysis=FALSE, seed=NULL){
     ### Define sources and simulation info
     which.sources <- as.data.table(t(sims.info[1:n.sources]))
     which.prop <- paste("H", n.sources, sep="")
     source.info <- cbind(which.sources, which.prop)
     which.sources <- unlist(which.sources)
     
     
     ### Define number of control objects, if necessary 
     if(length(n.control)==1){
          n.control <- rep(n.control, n.sources)
     }
     
     
     ### Sample spectra for control and trace objects 
     # Control objects
     if(!is.null(seed)) set.seed(seed)
     control.spectra <- NULL
     for(i in 1:n.sources){
          control.spectra <- cbind(control.spectra, 
                                   dat.ls[[which.sources[i]]][,(1:n.control[i])+1])
     }
     colnames(control.spectra) <- unlist(sapply(1:n.sources, function(x, n.control){
          paste("control.", x, '.', 1:n.control[x], sep="")
     }, n.control))
     # Test objects
     trace.spectra <- dat.ls[[which.sources[n.sources]]][,n.control[i]+1+(1:n.test)]
     if(!is.matrix(trace.spectra)){
          trace.spectra <- matrix(trace.spectra, ncol=1)
     }
     colnames(trace.spectra) <- paste("test.", 1:n.test.samps, sep="")
     # Combine the control and the trace spectra into one matrix
     spectra.mat <- cbind(control.spectra, trace.spectra)
     
     
     ### Get scores
     if(SVM.analysis==TRUE){
          ### Use common kernel for <method name> and SVM 
          all.scores <- apply(pairs.all, 1, function(x, spectra.mat){
               kernfun(spectra.mat[,x[1]], spectra.mat[,x[2]])
          }, spectra.mat)
          
          ### SVM 
          SVM.dat <- data.frame(Lab=rep(as.factor(which.sources), times=n.control), t(control.spectra))
          stop.crit = 1
          while(stop.crit <= 10){
             stop.crit = stop.crit + 1
             m.SVM <- ksvm(Lab~., data=SVM.dat, kernel=spectkern, scaled=FALSE, prob.model=TRUE, type="C-svc")
             tryCatch(p.SVM <- predict(m.SVM, t(trace.spectra), type='probabilities'), error = function(e) e)
             if(exists('p.SVM') == TRUE) stop.crit=11
          }
     } else{
          all.scores <- apply(pairs.all, 1, function(x, spectra.mat, xs.eval, max.lag, wav.num){
               kernfun(spectra.mat[,x[1]], spectra.mat[,x[2]])
          }, spectra.mat, xs.eval, max.lag, wav.num)
     }
     control.scores <- all.scores[unlist(ix.all[1:(n.sources+choose(n.sources, 2))])[control.ix]]
     
     
     ### Get Model Parameters via Adaptive MH
     # Define function for calculating log-posterior
     MCMC.logpost <- logpost.scores(control.scores, n.sources=n.sources, n.samps=n.control, PP.t.control) 
     param.samp.time <- system.time({
          param.samps <- MCMC(MCMC.logpost, init=theta0, adapt=TRUE, 
                              n=n.samps, acc.rate=0.234, 
                              scale=rep(1, 2*(n.sources+choose(n.sources,2))+1))$samples
     })
     param.samps.init = matrix(theta0, nrow=1)
     param.samps.tail <- param.samps[seq(0.75*nrow(param.samps), nrow(param.samps), 5),]
     
     ### Evaluate Bayes Factor Using Sampled Parameter Values 
     BF.eval.time <- system.time({
          BF.samps <- sapply(1:n.sources, BF.eval, param.samps.tail, all.scores, P.mat.all, n.sources, 
                             n.control, n.test.samps)
     })
     BF.num <- apply(BF.samps, 2, function(x) mean(exp(x[is.finite(x)]/100))) 
     BF <- sapply(BF.num, function(x, BF.num) x/sum(BF.num), BF.num)
     
     BF.samps.init <- matrix(sapply(1:n.sources, BF.eval, param.samps.init, all.scores, P.mat.all, n.sources, 
                        n.control, n.test.samps), nrow=1)
     BF.num.init <- apply(BF.samps.init, 2, function(x) mean(exp(x[is.finite(x)]/100))) 
     BF.init <- sapply(BF.num.init, function(x, BF.num.init) x/sum(BF.num.init), BF.num.init)
     
     
     
     ### Return relevant info
     if(SVM.analysis==TRUE){
          return(list(param.samps=param.samps, BF.samps=BF.samps, BF.samps.init=BF.samps.init, 
                      BF.num=BF.num, BF=BF, BF.num.init = BF.num.init, BF.init=BF.init,
                      n.sources=n.sources, param.samp.time=param.samp.time, BF.eval.time=BF.eval.time, 
                      p.SVM=p.SVM))
     } else{
          return(list(param.samps=param.samps, BF.samps=BF.samps, BF.samps.init=BF.samps.init, 
                      BF.num=BF.num, BF=BF, BF.num.init = BF.num.init, BF.init=BF.init,
                      n.sources=n.sources, param.samp.time=param.samp.time, BF.eval.time=BF.eval.time))
     }
}
################################################################################################# 

