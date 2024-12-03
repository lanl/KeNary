#################################################################################################
##### Look at results
#################################################################################################
### Load results
rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)


setwd('~/Documents/GitHub/kenary/')
source.grid <- as.matrix(data.table::fread(paste(getwd(), '/PaintSims/source_grid.csv', sep='')), nrow=100, ncol=10)
my.files <- list.files(getwd(), pattern='.rda')

# counter = 3
# tmp = list()
# # for(i in c(7,9:14)){
# for(i in c(15, 18:22)){
#   load(paste(getwd(), '/', my.files[i], sep=''))
#   tmp[[counter]] = FTIR.sims.results
#   counter = counter + 1
# }
# FTIR.sims.results = tmp
# save(FTIR.sims.results, file=paste(getwd(), '/FTIRSimsResults_8sources_3test_22Oct24.rda', sep=''))

my.files = my.files[grep(my.files, pattern="sources_3test")]
my.files = my.files[-c(6,8,9)]

# ### Load a test case
# which.file =1
# load(paste(getwd(), '/', my.files[which.file], sep=''))                                                        
# # NOTES: 
# # opens as FTIR.sims.results
# # list of length 10 
# # 100 elements in each list element
# 
# n.source = as.numeric(substr(paste(getwd(), '/', my.files[which.file], sep=''), 59, 59))
# KBSC.rate <- SVM.rate <- list()
# for(i in 3:10){
#  SVM.rate[[i]] <- matrix(0, nrow=length(FTIR.sims.results[[i]]), ncol=3)
#  KBSC.rate[[i]] <- sum(unlist(lapply(FTIR.sims.results[[i]], function(x) which.max(x$BF)==x$n.sources)))/length(FTIR.sims.results[[i]])
#  for(j in 1:length(FTIR.sims.results[[i]])){
#   tmp.col <- apply(FTIR.sims.results[[i]][[j]]$p.SVM, 1, which.max)
#   SVM.rate[[i]][j,] <- names(sapply(1:3, function(x, sims.results, tmp.col) 
#                              sims.results[x, tmp.col[x]], FTIR.sims.results[[i]][[j]]$p.SVM, tmp.col)) == source.grid[j, n.source]
#  }
# }
# 
# # Look at KBSC Success Rate
# unlist(KBSC.rate)
# 
# # Look at SVM Success Rate 
# # Get per-sim success 
# sim.rates <- lapply(SVM.rate[3:length(SVM.rate)], rowMeans)
# # all or nothing 
# all.or.nothing <- lapply(sim.rates, function(x) ifelse(x==1, 1, 0))
# unlist(lapply(all.or.nothing, mean))
# # per-object 
# per.object <- lapply(sim.rates, mean)
# unlist(per.object)
# # voting 
# voting <- lapply(sim.rates, function(x) ifelse(x>=2/3, 1, 0))
# unlist(lapply(voting, mean))
# 
# # Look at Rhat 
# load(paste(getwd(), '/', my.files[2], sep=''))                                                        
# 
# rhat <- list()
# for(i in 3:length(FTIR.sims.results)){
#  rhat[[i]] <- lapply(FTIR.sims.results[[i]], function(x){
#   tmp <- rep(0, ncol(x$param.samps))
#   for(j in 1:ncol(x$param.samps)){
#    tmp[j] <- (rstan::Rhat(x$param.samps[seq(7000,length(x$param.samps[,j]), 50),j]))
#   }
#   return(tmp)
#  })
#  rhat[[i]] <- unlist(lapply(rhat[[i]], median))
# }
# lapply(rhat, median, na.rm=TRUE)

#################################################################################################
##### Presenting Results 
#################################################################################################

# 2 sources, 3 control, 3 test 
facet_names <- list(
  '1'=bquote(theta[11]),
  '2'=bquote(theta[22]),
  '3'=bquote(theta[33]),
  '4'=bquote(theta[12]),
  '5'=bquote(theta[13]),
  '6'=bquote(theta[23]),
  '7'=bquote(sigma[11]),
  '8'=bquote(sigma[22]),
  '9'=bquote(sigma[33]),
  '10'=bquote(sigma[12]),
  '11'=bquote(sigma[13]),
  '12'=bquote(sigma[23])
)
facet_labeller <- function(variable,value){
  return(facet_names[value])
}

which.mix = 78

load(paste(getwd(), '/', my.files[1], sep=''))
my.seq <- seq(0.70*nrow(FTIR.sims.results[[3]][[which.mix]]$param.samps), nrow(FTIR.sims.results[[3]][[which.mix]]$param.samps), by=10)
for(i in c(1,2)){
  
  load(paste(getwd(), '/', my.files[i], sep=''))
  this.min=nrow(FTIR.sims.results[[3]][[which.mix]]$param.samps)
  n.cols = ncol(FTIR.sims.results[[3]][[which.mix]]$param.samps) - 1
  tmp.df = reshape2::melt(abs(FTIR.sims.results[[3]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
  if(i==1){
    par.vect = c(58.70, 68.27, 85.01, 2.61, 5.26, 0.27)
    # par.vect = c(71.50, 74.30, 87.31, 3.26, 3.09, 1.14) #57
    mean.vect = colMeans(abs(FTIR.sims.results[[3]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
    median.vect = apply(abs(FTIR.sims.results[[3]][[which.mix]]$param.samps[my.seq, 1:n.cols]), 2, median)
    par.df = data.frame(Var1=par.vect, Var2=1:6, Var3=mean.vect, Var4=median.vect)
  }else if(i==2){
    par.vect = c(58.70, 68.27, 63.52, 85.01, 79.75, 84.72, 2.61, 5.26, 1.84, 0.27, 0.49, 0.32)
    #par.vect = c(71.50, 74.30, 69.78, 87.31, 93.50, 91.06, 3.26, 3.09, 3.09, 1.14, 0.38, 0.45) #57
    mean.vect = colMeans(abs(FTIR.sims.results[[3]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
    median.vect = apply(abs(FTIR.sims.results[[3]][[which.mix]]$param.samps[my.seq, 1:n.cols]), 2, median)
    par.df = data.frame(Var1=par.vect, Var2=1:12, Var3=mean.vect, Var4=median.vect)
  }
  
  MH3 = ggplot(tmp.df, aes(x=Var1, y=value)) + 
    ylab("") +
    xlab("") +
    geom_line() + 
    geom_hline(data=par.df, aes(yintercept=Var1, col="Empirical Value")) +
    geom_hline(data=par.df, aes(yintercept=Var3, col="Chain Mean")) +
    geom_hline(data=par.df, aes(yintercept=Var4, col="Chain Median")) +
    facet_wrap(vars(Var2), scales='free_y', labeller=facet_labeller, nrow=2) + 
    #guides(color=guide_legend(title="Statistics")) +
    theme(legend.position="", 
          legend.text=element_text(size=15), 
          legend.title=element_text(size=20), 
          axis.text.x=element_text(size=15), 
          axis.text.y=element_text(size=15), 
          strip.text=element_text(size=17), 
          title=element_text(size=20)) + 
    labs(title=paste("Adaptive MH Samples for Combination ", which.mix, ": ", i+1, " sources, 3 control objects", sep="")) + 
    scale_y_continuous(n.breaks=5)
  
  rhat.all <- rhat.mean <- list()
  rhat.all <- lapply(FTIR.sims.results[[3]], function(x){
    tmp <- rep(0, ncol(x$param.samps))
    for(j in 1:ncol(x$param.samps)){
      tmp[j] <- (rstan::Rhat(x$param.samps[my.seq,j]))
    }
    return(tmp)
  })
  
  rhat.mean <- unlist(lapply(rhat.all, median))
  print(paste('Rhat:', rhat.all[[which.mix]]))
  print(rhat.mean[which.mix])


  tmp.df = reshape2::melt(abs(FTIR.sims.results[[8]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
  if(i==1){
    par.vect = c(65.83, 66.90, 84.83, 5.93, 4.78, 0.36)
    # par.vect = c(66.54, 77.71, 98.52, 7.52, 4.33, 0.60) #52
    mean.vect = colMeans(abs(FTIR.sims.results[[10]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
    median.vect = apply(abs(FTIR.sims.results[[10]][[which.mix]]$param.samps[my.seq, 1:n.cols]), 2, median)
    par.df = data.frame(Var1=par.vect, Var2=1:6, Var3=mean.vect, Var4=median.vect)
  }else{
    par.vect = c(65.83, 66.90, 70.47, 84.83, 81.14, 85.07, 5.93, 4.78, 6.41, 0.36, 2.35, 0.93)
    # par.vect = c(68.15 76.33, 68.55, 87.95, 93.41, 91.06, 4.73, 5.71, 4.81, 2.02, 0.40, 0.57) #10 control objects
    mean.vect = colMeans(abs(FTIR.sims.results[[10]][[which.mix]]$param.samps[my.seq, 1:n.cols]))
    median.vect = apply(abs(FTIR.sims.results[[10]][[which.mix]]$param.samps[my.seq, 1:n.cols]), 2, median)
    par.df = data.frame(Var1=par.vect, Var2=1:12, Var3=mean.vect, Var4=median.vect)
  }
  MH8 = ggplot(tmp.df, aes(x=Var1, y=value)) +
    xlab("Adpative MH Iteration") +
    ylab("") +
    geom_line() +
    geom_hline(data=par.df, aes(yintercept=Var1, col="Empirical Value")) +
    geom_hline(data=par.df, aes(yintercept=Var3, col="Chain Mean")) +
    geom_hline(data=par.df, aes(yintercept=Var4, col="Chain Median")) +
    facet_wrap(vars(Var2), scales='free_y', labeller=facet_labeller, nrow=2) +
    guides(color=guide_legend(title="Statistics")) +
    theme(legend.position="bottom",
          legend.text=element_text(size=15),
          legend.title=element_text(size=20),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          strip.text=element_text(size=17),
          title=element_text(size=20)) +
    labs(title=paste("Adaptive MH Samples for Combination ", which.mix, ": ", i+1, " sources, 8 control objects", sep="")) + 
    scale_y_continuous(n.breaks=5)

  rhat.all <- rhat.mean <- list()
  rhat.all <- lapply(FTIR.sims.results[[8]], function(x){
    tmp <- rep(0, ncol(x$param.samps))
    for(j in 1:ncol(x$param.samps)){
      tmp[j] <- (rstan::Rhat(x$param.samps[my.seq,j]))
    }
    return(tmp)
  })

  # ess.all <- lapply(FTIR.sims.results[[10]], function(x){
  #   tmp <- rep(0, ncol(x$param.samps))
  #   for(j in 1:ncol(x$param.samps)){
  #     tmp[j] <- (rstan::ess_bulk(x$param.samps[my.seq,j]))
  #   }
  #   return(tmp)
  # })

  rhat.mean <- unlist(lapply(rhat.all, median))
  print(paste('Rhat:', rhat.all[[which.mix]]))
  print(rhat.mean[which.mix])

  # ess.mean <- unlist(lapply(ess.all ,median))
  # print(paste('ESS:', ess.all[[which.mix]]))
  #print(ess.mean[which.mix])

  grid.arrange(MH3, MH8)
  
  # grid.arrange(MH3)
}



#####################################
#### 2 Class: KeNary vs. SVM 
#####################################
### Get KBSC success rate and SVM success rate for all numbers of classes
KBSC.rate = SVM.rate = sim.rates = all.or.nothing = per.object = voting = list()
for(which.file in 1:length(my.files)){
  load(paste(getwd(), '/', my.files[which.file], sep=''))                                                        
  n.source = as.numeric(substr(paste(getwd(), '/', my.files[which.file], sep=''), 59, 59))
  
  KBSC.rate[[which.file]] <- SVM.rate[[which.file]] <- list()
  for(i in 3:length(FTIR.sims.results)){
    SVM.rate[[which.file]][[i]] <- matrix(0, nrow=length(FTIR.sims.results[[i]]), ncol=3)
    KBSC.rate[[which.file]][[i]] <- sum(unlist(lapply(FTIR.sims.results[[i]], function(x) which.max(x$BF)==x$n.sources)))/length(FTIR.sims.results[[i]])
    for(j in 1:length(FTIR.sims.results[[i]])){
      tmp.col <- apply(FTIR.sims.results[[i]][[j]]$p.SVM, 1, which.max)
      SVM.rate[[which.file]][[i]][j,] <- names(sapply(1:3, function(x, sims.results, tmp.col) 
        sims.results[x, tmp.col[x]], FTIR.sims.results[[i]][[j]]$p.SVM, tmp.col)) == source.grid[j, n.source]
    }
  }
  
  # Look at SVM Success Rate 
  # Get per-sim success 
  sim.rates[[which.file]] <- lapply(SVM.rate[[which.file]][3:length(FTIR.sims.results)], function(x) rowMeans(x))
  
  # all or nothing 
  all.or.nothing[[which.file]] <- lapply(sim.rates[[which.file]] , function(x) ifelse(x==1, 1, 0))
  all.or.nothing[[which.file]] = lapply(all.or.nothing[[which.file]], mean)
  
  # per-object 
  per.object[[which.file]] <- lapply(sim.rates[[which.file]], mean)
  
  #voting 
  voting[[which.file]] <- lapply(sim.rates[[which.file]], function(x) ifelse(x>=2/3, 1, 0))
  voting[[which.file]] = lapply(voting[[which.file]], mean)
}


### Look at KBSCS for 2 classes as number of control objects varies 3:10
KBSC.df <- do.call('rbind', lapply(KBSC.rate, unlist))
#KBSC.df <- data.frame(KBSC=KBSC.df)
colnames(KBSC.df) = paste(3:10, 'control obs')
rownames(KBSC.df) = paste((1:length(my.files))+1, 'sources')
KBSC.df = reshape::melt(KBSC.df)
KBSC.df = data.frame(KBSC.df, method=rep('KeNary', length(my.files)*8))

### Look at SVM results for 2 classes as number of control objects varies 3:10
# all or nothing
all.or.nothing.df <- do.call('rbind', lapply(all.or.nothing, unlist))
#all.or.nothing.df <- as.data.frame(AorN=all.or.nothing.df)
colnames(all.or.nothing.df) = paste(3:10, 'control obs')
rownames(all.or.nothing.df) = paste((1:length(my.files))+1, 'sources')
all.or.nothing.df = reshape::melt(all.or.nothing.df, id.vars=1:2)
all.or.nothing.df = data.frame(all.or.nothing.df, method=rep('SVM (all or nothing)', length(my.files)*8))
# per object
per.object.df <- do.call('rbind', lapply(per.object, unlist))
#per.object.df <- data.frame(PO=per.object.df)
colnames(per.object.df) = paste(3:10, 'control obs')
rownames(per.object.df) = paste((1:length(my.files))+1, 'sources')
per.object.df = reshape::melt(per.object.df)
per.object.df = data.frame(per.object.df, method=rep('SVM (per object)', length(my.files)*8))
# voting 
voting.df <- do.call('rbind', lapply(voting, unlist))
#voting.df <- data.frame(VOTE=voting.df)
colnames(voting.df) = paste(3:10, 'control obs')
rownames(voting.df) = paste((1:length(my.files))+1, 'sources')
voting.df = reshape::melt(voting.df)
voting.df = data.frame(voting.df, method=rep('SVM (voting)', length(my.files)*8))

### Combine all results into a single data frame
all.results.df <- rbind(KBSC.df, all.or.nothing.df, per.object.df, voting.df)
all.results.df <- data.frame(all.results.df %>% dplyr::group_by(method, X1) %>% dplyr::mutate(Avg_Success=mean(value)))
levels(all.results.df$X2) = c(3:10)

# ggplot(all.results.df, aes(x=X1, y=value, color=method)) + 
#   geom_boxplot(alpha=0.5, aes(fill=method), position=position_dodge(1.)) + 
#   #geom_point(aes(alpha=X2, group=method), position=position_jitterdodge(jitter.width=0.5), size=2.5) +
#   scale_alpha_manual(values=seq(0.25, 1, length.out=8)) +
#   ggthemes::scale_color_colorblind(guide="none") +
#   ggthemes::scale_fill_colorblind() +
#   theme(legend.position='bottom', legend.box='vertical', 
#         axis.text.x=element_text(size=20), axis.title.x=element_text(size=25), 
#         axis.text.y=element_text(size=20), axis.title.y=element_text(size=25),
#         legend.text=element_text(size=20), legend.title=element_text(size=25), title=element_text(size=28)) + 
#   labs(title="Performance of KeNary vs. SVM", x="Number of Sources", y="Percent Correct Classification") + 
#   guides(fill=guide_legend(title="Classification Method")) +
#   ggrepel::geom_label_repel(aes(label=X2, group=method), 
#                             position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0.0001, dodge.width = 1.), 
#                             max.time=10, direction='both', label.padding=0.15, label.size=0.4, size=8)

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
  #facet_wrap(~X1, ncol=1) 


design <- c(
  "
#AABBCC#
DDEEFFGG
"
)

p + ggh4x::facet_manual(~X1, design = design)

  #ggrepel::geom_label_repel(aes(label=X2, group=method), 
  #                          position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0.0001, dodge.width = 1.), 
  #                          max.time=10, direction='both', label.padding=0.15, label.size=0.4, size=8)


# ggplot(all.results.df, aes(x=X1, y=value, fill=method, alpha=X2)) +
#   geom_col(position='dodge') +
#   scale_alpha_manual(values=seq(0.5, 1, length.out=8)) +
#   theme(legend.position='bottom', legend.box='vertical', text=element_text(size=15)) +
#   ggthemes::scale_fill_colorblind(guide="none") +
#   #coord_cartesian(ylim = c(0.85, 1)) +
#   geom_line(aes(x=X1, y=Avg_Success, col=method, group=method)) +
#   ggthemes::scale_colour_colorblind(guide='none') + facet_wrap(~method) +
#   labs(title="Performance of KeNary vs. SVM", x="Number of Sources", y="Percent Correct Classification") +
#   guides(alpha=guide_legend(title="Number of Control Objects"))


levels(all.results.df$X2) = c('3', '4', '5', '6', '7', '8', '9', '10')
ggplot(subset(all.results.df, !(X2%in%c('9', '10'))), 
       aes(x=X2, y=value, fill=method, shape=method, color=method)) + 
  geom_point(position=position_dodge(width=0.35), size=4, stroke=1.) +
  scale_shape_manual(values=c(23, 21,21,21)) + 
  scale_color_manual(values=c('red3',1,1,1)) +
  theme(legend.position='bottom', legend.box='vertical', 
        axis.text.x=element_text(size=18), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=18), axis.title.y=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
        strip.text=element_text(size=18)) + 
  coord_cartesian(ylim=c(0.60, 1.015)) + 
  #geom_line(aes(group=method), position=position_jitterdodge()) + 
  ggthemes::scale_fill_colorblind() +
  facet_wrap(~X1, ncol=1, strip.position="top") + 
  #scale_x_discrete(labels=c('3', '4', '5', '6', '7', '8'), breaks=3:8) +
  labs(title="Performance of KeNary vs. SVM with 3 to 8 Control Objects", x="Number of Control Objects", y="Percent Correct Classification") + 
  guides(fill=guide_legend(title="Classification Method"), shape=guide_legend(title="Classification Method"), color=guide_legend(title="Classification Method"))


load(paste(getwd(), '/', my.files[1], sep=''))
tmp=list()
for(i in 3:10){
  tmp.df = exp(FTIR.sims.results[[i]][[which.mix]]$BF.samps/10000)
  tmp.rowsum = rowSums(tmp.df)
  
  tmp.df[,1] = tmp.df[,1]/tmp.rowsum
  tmp.df[,2] = tmp.df[,2]/tmp.rowsum
  tmp[[i]] = reshape2::melt(tmp.df)[,-1]
  tmp[[i]] = cbind(tmp[[i]], control=rep(i, nrow(tmp[[i]])))
}
tmp.df = do.call('rbind', tmp)
ggplot(subset(tmp.df, control%in%c(3,8)), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
  stat_density(alpha=0.5) + 
  scale_fill_manual(values=c("#F0E442", '#0072B2')) + 
  scale_color_manual(values=c("#F0E442", '#0072B2'), guide="none") + 
  guides(fill=guide_legend(title="Probability of Class... "), col=guide_legend(title="Probability of Class... ")) + 
  labs(title=paste("Distribution of Posterior Probability of Class for Combination ", which.mix, ' (2 Sources)')) + 
  xlab("") +
  theme(legend.position='bottom', legend.box='vertical', 
        axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
        strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  facet_wrap(~control, scales="free", nrow=1) + 
  scale_x_continuous(n.breaks=3) 

source(paste(getwd(), '/AdaptiveMH_DarwinPaintFuns.R', sep=''))
BF.samps.init <- matrix(sapply(1:2, BF.eval, FTIR.sims.results[[3]][[which.mix]]$param.samps[1:1000,], all.scores, P.mat.all, 2, 
                               3, 3), nrow=1)
BF.num.init <- apply(BF.samps.init, 2, function(x) mean(exp(x[is.finite(x)]/100))) 
BF.init <- sapply(BF.num.init, function(x, BF.num.init) x/sum(BF.num.init), BF.num.init)

tmp=list()
for(i in 3:10){
  tmp.df = exp(FTIR.sims.results[[3]][[which.mix]]$BF.samps)
  tmp.rowsum = rowSums(tmp.df)
  
  tmp.df[,1] = tmp.df[,1]/tmp.rowsum
  tmp.df[,2] = tmp.df[,2]/tmp.rowsum
  tmp[[i]] = reshape2::melt(tmp.df)[,-1]
  tmp[[i]] = cbind(tmp[[i]], control=rep(i, nrow(tmp[[i]])))
}
tmp.df = do.call('rbind', tmp)
ggplot(subset(tmp.df, control%in%c(3,8)), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
  stat_density(alpha=0.5) + 
  scale_fill_manual(values=c("#F0E442", '#0072B2')) + 
  scale_color_manual(values=c("#F0E442", '#0072B2'), guide="none") + 
  guides(fill=guide_legend(title="Probability of Class... "), col=guide_legend(title="Probability of Class... ")) + 
  labs(title=paste("Distribution of Prior Probability of Class for Combination ", which.mix, ' (2 Sources)')) + 
  xlab("") +
  theme(legend.position='bottom', legend.box='vertical', 
        axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
        strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  facet_wrap(~control, scales="free", nrow=1) + 
  scale_x_continuous(n.breaks=3) 

# prior plots
theta0 = c(rep(50,choose(s,2)+s), rep(5, choose(s,2)+s), 0.1)
theta0.mat=matrix(runif(25*7, -0.5, 0.5), nrow=25) + matrix(theta0, nrow=25, ncol=7, byrow=T)
dbeta(theta0[length(theta0)], 0.5, 0.5, log=TRUE) + 
  sum(dlomax(theta0[1:(n.sources+choose(n.sources,2))], shape=6, scale=1)) +
  sum(dhalfcauchy(theta0[n.sources+choose(n.sources,2)+1:(n.sources+choose(n.sources,2))], 1e-5, log=TRUE))

my.prior.fun = function(par){
  dbeta(par[length(par)], 0.5, 0.5, log=TRUE) + 
  sum(dlomax(par[1:(n.sources+choose(n.sources,2))], shape=6, scale=1)) +
  sum(dhalfcauchy(par[n.sources+choose(n.sources,2)+1:(n.sources+choose(n.sources,2))], 1e-5, log=TRUE))
}
apply(theta0.mat, 1, my.prior.fun)


load(paste(getwd(), '/', my.files[2], sep=''))
for(i in 3:10){
  tmp.df = exp(FTIR.sims.results[[i]][[which.mix]]$BF.samps/25000)
  tmp.rowsum = rowSums(tmp.df)
  
  tmp.df[,1] = tmp.df[,1]/tmp.rowsum
  tmp.df[,2] = tmp.df[,2]/tmp.rowsum
  tmp.df[,3] = tmp.df[,3]/tmp.rowsum
  #tmp.df[,4] = tmp.df[,4]/tmp.rowsum
  #tmp.df[,5] = tmp.df[,5]/tmp.rowsum
  tmp[[i]] = reshape2::melt(tmp.df)[,-1]
  tmp[[i]] = cbind(tmp[[i]], control=rep(i, nrow(tmp[[i]])))
}
tmp.df = do.call('rbind', tmp)

#### Make sure above sums to 1!!!! 

ggplot(subset(tmp.df, control%in%c(3,8)), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
  stat_density(alpha=0.5) + 
  scale_fill_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")[1:3]) + 
  scale_color_manual(values=c("#F0E442", "#0072B2", "#CC79A7", "#009E73", "#D55E00")[1:3], guide="none") + 
  guides(fill=guide_legend(title="Probability of Class... ")) + 
  labs(title=paste("Distribution of Probability of Class for Combination ", which.mix, ' (3 Sources)', sep="")) + 
  xlab("") +
  theme(legend.position='bottom', legend.box='vertical', 
        axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
        strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  facet_wrap(~control, scales="free", nrow=1) + 
  scale_x_continuous(n.breaks=3)






load(paste(getwd(), '/', my.files[4], sep=''))
for(i in 3:10){
  tmp.df = exp(FTIR.sims.results[[i]][[which.mix]]$BF.samps/150000)
  tmp.rowsum = rowSums(tmp.df)
  
  tmp.df[,1] = tmp.df[,1]/tmp.rowsum
  tmp.df[,2] = tmp.df[,2]/tmp.rowsum
  tmp.df[,3] = tmp.df[,3]/tmp.rowsum
  tmp.df[,4] = tmp.df[,4]/tmp.rowsum
  tmp.df[,5] = tmp.df[,5]/tmp.rowsum
  tmp[[i]] = reshape2::melt(tmp.df)[,-1]
  tmp[[i]] = cbind(tmp[[i]], control=rep(i, nrow(tmp[[i]])))
}
tmp.df = do.call('rbind', tmp)


ggplot(subset(tmp.df, control%in%c(3,8)), aes(x=value, col=as.factor(Var2), group=as.factor(Var2), fill=as.factor(Var2))) + 
  stat_density(alpha=0.5) + 
  scale_fill_manual(values=c("#F0E442", '#0072B2', "#CC79A7", "#009E73", "#D55E00")) + 
  scale_color_manual(values=c("#F0E442", "#0072B2", "#CC79A7", "#009E73", "#D55E00"), guide="none") + 
  guides(fill=guide_legend(title="Probability of Class... ")) + 
  labs(title=paste("Distribution of Probability of Class for Combination ", which.mix, ' (5 Sources)', sep="")) + 
  xlab("") +
  theme(legend.position='bottom', legend.box='vertical', 
        axis.text.x=element_text(size=13), axis.title.x=element_text(size=20), axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=13), axis.title.y=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20), title=element_text(size=21), 
        strip.text=element_text(size=18), plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  facet_wrap(~control, scales="free", nrow=1) + 
  scale_x_continuous(n.breaks=3)





