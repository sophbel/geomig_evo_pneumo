library(data.table)
library(Rcpp)
# 
# ######################################################################################################
# ######MUNICIPALITY LEVEL##########
# #######################################################################################################
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
load("./modelinput_data/pairwise_geodist.town.RData") # # [21_05_21_PariwiseDistance.R] ## Pairwise distance between towns in South Africa
load("./modelinput_data/pop_municipality.2017LS.RData")
load("./modelinput_data/pop_2019.RData")
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/tmpbase.RData")
load("./modelinput_data/cdr.mat.town.one.RData")# # [Mobility_ManyMonthsSA.R] ## Probability of movement between each province normalized to carriage rates for each province
cdr.mat.town<-cdr.mat.town.one
tn <- rownames(cdr.mat.town)
# tmpbase<-tmpbase[tn,tn]
pairwise_geodist.town<-pairwise_geodist.town[tn,tn]
pop2019.town <-pop2019.town[tn]
nlocs=dim(pairwise_geodist.town)[1]
distmin <- c(0,10,200,500,1000)
distmax <- c(10,200,500,1000,2000)
genvec <- c(1,10)
# genvec <- c(1,6)
load("./ModelProjections/data/TranMatArray.234x234_adj.RData")
TranMatArray.1<-TranMatArray.234x234
# 

# ######BOOTSTRAP CONFIDENCE INTERVALS##########
# TranMA.tmp <- array(data=NA, dim=c(nlocs,nlocs,length(genvec)))
# iters=1000
# pair.rand <- matrix(1/234,nrow=234,ncol=234)
# dist.prop.mat.list <-list()
# for (m in 1:iters){
#   dist.prop.mat <- matrix(nrow=length(genvec), ncol=length(distmax))
#   dist.prop.mat.rand <-vector(mode="numeric",length=length(distmax))
#   for (j in 1:length(genvec)){
#     prop.mat<-matrix(nrow=length(distmax),ncol=nlocs)
#     prop.mat.rand<-matrix(nrow=length(distmax),ncol=nlocs)
#     for (i in 1:nlocs){
#       # startSamp <- sample(234,1,replace=T,prob=mrcaVec)
#       startSamp <- sample(234,1,replace=T,prob=pop2019.town)
#       
#       # startSamp <- sample(234,1,replace=T)
#     for (k in 1:length(distmax)){
#       prop.mat[k,i] <- sum(TranMatArray.1[startSamp,,genvec[j]][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
#       prop.mat.rand[k,i]<-sum(pair.rand[startSamp,][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
#       }
#       # print(i)
#     }
#     dist.prop.mat[j,] <- apply(prop.mat,1,function(x) mean(x))
#     dist.prop.mat.rand  <- apply(prop.mat.rand,1,function(x) mean(x))
#     # dist.prop.mat[j,] <- apply(prop.mat,1,function(x) weighted.mean(x,pop2019.town))
#     # dist.prop.mat.rand  <- apply(prop.mat.rand,1,function(x) weighted.mean(x,pop2019.town))
#   }
#   dist.prop.mat <- data.table(dist.prop.mat)
#   dist.prop.mat.list[[m]] <- dist.prop.mat
#   print(m)
# }
# 
# 
# 
# ### Find quantiles across bootstraps
# dist.prop.mat2.tmp <- rbindlist(dist.prop.mat.list)
# colnames(dist.prop.mat2.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
# dist.prop.mat2.tmp$gens <- rep(as.character(genvec),iters)
# dist.prop.mat.rand.tmp<-data.table(t(data.table(as.numeric(dist.prop.mat.rand))))
# colnames(dist.prop.mat.rand.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
# dist.prop.mat.rand.tmp$gens <- "Random"
# 
# dist.prop.mat2<-rbind(dist.prop.mat2.tmp,dist.prop.mat.rand.tmp)
# dist.prop.mat2 <- melt(dist.prop.mat2)
# colnames(dist.prop.mat2) <- c("gens","distances","value")
# dist.prop.mat2$distances <- factor(dist.prop.mat2$distances,levels=c("0-10km","10-200km","200-500km","500-1000km",">1000km"),labels=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"))
# genvec2<-c(1,10,50,"Random")
# # genvec2<-c(1,6,"Random")
# 
# ## Quantiles ### DOTS AND REDUCE THE NUMBER OF GENERATIONS
# dists <- c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km")
# dist.quants.list<-list()
# for (i in 1:length(dists)) {
#   dist.quants.tmp <- matrix(nrow=length(genvec2),ncol=6)
#   for (j in 1:length(genvec2)){
#     tmp <- subset(dist.prop.mat2, dist.prop.mat2$distances == dists[i] & dist.prop.mat2$gens==genvec2[j])
#     dist.quants.tmp[j,1:3] <- as.numeric(quantile(tmp$value, probs=c(0.025,0.5,0.975)))
#     dist.quants.tmp[j,4] <- as.numeric(mean(tmp$value))
#     dist.quants.tmp[j,5] <- dists[i]
#     dist.quants.tmp[j,6] <- genvec2[j]
#   }
#   dist.quants.tmp <- data.table(dist.quants.tmp)
#   dist.quants.list[[i]] <- dist.quants.tmp
# }
# dist.quants <- rbindlist(dist.quants.list)
# colnames(dist.quants) <- c("lowerCI","median","upperCI","mean","distances","gens")
# dist.quants$lowerCI <- as.numeric(dist.quants$lowerCI)
# dist.quants$median <- as.numeric(dist.quants$median)
# dist.quants$mean <- as.numeric(dist.quants$mean)
# 
# dist.quants$upperCI <- as.numeric(dist.quants$upperCI)
# 
# # dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==6 | dist.quants$gens=="Random")
# dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==10 | dist.quants$gens==50 |dist.quants$gens=="Random")
# dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==10 | dist.quants$gens=="Random")
# 
# dist.quants.tmp3$lowerCI[dist.quants.tmp3$gens=="Random"] <- NA
# dist.quants.tmp3$upperCI[dist.quants.tmp3$gens=="Random"] <- NA
# 
# ggplot(data=dist.quants.tmp3) +
#   geom_point(data=dist.quants.tmp3,aes(x=distances,y=median,shape=gens),size=3,position=position_dodge(width=0.35)) +
#   geom_errorbar(aes(ymin=lowerCI,ymax=upperCI,shape=gens,x=distances),alpha=0.6,width=0.2,position=position_dodge(width=0.35))+
#   xlab("Distance \n(from starting municipality)")+
#   ylab("Proportion")+
#   theme(panel.grid.minor  = element_blank(),)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90,size=18), axis.text.y=element_text(size=18),
#         axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
#         legend.text = element_text(size=12),legend.title = element_text(size=12),
#         legend.position = c(.8,.8),legend.background = element_rect(color="black", size=0.5, linetype="solid"))+
#   # scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,6,"Random"), labels=c("1 Generation","1 Year","Given random \nmobility")) +
#   # scale_shape_manual(name="Time",values=c(19,2,8,7),limits=c(1,10,50,"Random"), labels=c("1 Gen","10 Gens","50 Gens","Given random \nmobility")) +
#   scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,10,"Random"), labels=c("1 Gen","10 Gens","Given random \nmobility")) +
#   scale_x_discrete(limits=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"),labels=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"))+
#   guides(fill="none") +
#   ylim(0,1)
# 
# # ggsave("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/propPerDist.sampUncertainty.pdf",width=4,height=5.5)
# 
# 








####Mean Across any location --- parameters for confidence intervals 
extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999
### all chains
# setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/RunModel/outputs/munic/")
load("./MCMC_model/outputs/mechan_model/ans.munic1.20000.08_adj.RData")
chain1<-ans.munic
# chain1<-mcmc(chain1,start=3000,end=4000)
chain1<-chain1[3000:4000,]
load("./MCMC_model/outputs/mechan_model/ans.munic2.20000.08_adj.RData")
chain2<-ans.munic
# chain2<-mcmc(chain2,start=3000,end=4000)
chain2<-chain2[3000:4000,]
load("./MCMC_model/outputs/mechan_model/ans.munic3.20000.08_adj.RData")
chain3<-ans.munic
# chain3<-mcmc(chain3,start=3000,end=4000)
chain3<-chain3[3000:4000,]
ans.munic.chains<-rbind(chain1,chain2,chain3)

ans<-ans.munic.chains
posteriors<-as.matrix(ans)
nsim<-1000
pars<-matrix(nrow=nsim,ncol=ncol(ans))
for(i in 1:nsim){
  simno<-sample(nrow(posteriors),1)
  pars[i,]<-posteriors[simno,]
}
sourceCpp("./MCMC_model/MatrixMultiplication.cpp")

Eigencpp=TRUE
maxGen=12
tn <- rownames(cdr.mat.town)
pairwise_geodist.town<-pairwise_geodist.town[tn,tn]
pop2019.town <-pop2019.town[tn]
nlocs=dim(pairwise_geodist.town)[1]
distmin <- c(0,10,200,500,1000)
distmax <- c(10,200,500,1000,2000)
genvec <- c(1,10)


## Load in the sampled parameters to bootstrap
dist.daily_probs.list<-dist.prop.mat.list <-dist.propSingle.list<-list()
dist.prop.mat.rand <-dist.propFB<-matrix(nrow=234,ncol=length(distmax))
pair.rand <- matrix(1/234,nrow=234,ncol=234)

### ###Load new mobility mat
min.range<-(-0.04)
max.range<-0.6
tranmat.cdr<-cdr.mat.town[tn,tn] #cdr.mat.town # My 9X9 matrix
### Many ways of altering parameter value
  for (m in 1:nrow(pars)){
    # for (m in 1:10){ 
    # fit.par<-apply(ans.munic.chains,2,median)
    fit.par<-pars[m,]
  # fit.par<-summary(ans)$statistics[,1]
  # fit.par[1]
  # nInfecLoc<-rep(1,9)
  # nInfecLoc[1:8]<-exp(fit.par[2:9])
  # nInfecLoc<-nInfecLoc/sum(nInfecLoc)
  # probByloc<-no.detect/nInfecLoc
  # 
  
  nloc=nlocs=9
  tmp.pHome<-pop2019.town

  extTranMatDat.tmp$pars$homeSus<-fit.par[1]
  
  ### Create mobility matrix for susceptible individuals
  tmpbase_pre<-cdr.mat.town
  nloc.munic<-nrow(cdr.mat.town)
  tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
  tmppar <- min.range+tmppar1*(max.range-min.range)
  tmpdiag<-diag(tmpbase_pre)-tmppar
  tmpdiag[which(tmpdiag>0.99999)]<-0.99999
  diag(tmpbase_pre)<-0
  tmpbase_pre1<-sweep(tmpbase_pre,1,rowSums(tmpbase_pre),"/")
  tmpbase_pre2<-sweep(tmpbase_pre1,1,(1-tmpdiag)/(1-diag(tmpbase_pre1)),"*")
  diag(tmpbase_pre2)<-tmpdiag
  mean(tmpdiag)
  quantile(tmpdiag,probs=c(0.025,0.5,0.975))
  
  #### adjust so it is mobility across the infectious period
  timeWindow<-35
  probStay<-1-(diag(tmpbase_pre2))^timeWindow
  tmp<-tmpbase_pre2
  diag(tmp)<-0
  tmp1<-sweep(tmp,1,rowSums(tmp),"/")
  tmp2<-sweep(tmp1,1,(1-probStay)/(1-diag(tmp1)),"*")
  diag(tmp2)<-probStay
  tmpbase<-tmp2
  quantile(diag(tmpbase),probs=c(0.025,0.5,0.975))
  
  move3<-tcrossprod(tmpbase,tmpbase)
  move4<-sweep(move3,2,tmp.pHome,"*")
  TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")
  
  mean(diag(TranMat.tmp2))
  weighted.mean(diag(TranMat.tmp2),pop2019.town)
  
  probByGenA<-probByGenB<-rep(1,maxGen)
  gensB<-gensA<-(1:maxGen)
  maxgen.tmp<-max(c(gensA,gensB))
  TranMatArray.1<-array(NA,c(nloc.munic,nloc.munic,maxgen.tmp))
  TranMatArray.1[,,1]<-TranMat.tmp2
  if(Eigencpp){
    for (j in 2:maxgen.tmp){TranMatArray.1[,,j]<-eigenMapMatMult(TranMatArray.1[,,j-1], TranMat.tmp2)
    }
  }else{
    for (j in 2:maxgen.tmp){TranMatArray.1[,,j]<-TranMatArray.1[,,j-1]%*% TranMat.tmp2}
  }
  
  locs=234
prop.mat.gens <- matrix(nrow=length(genvec),ncol=length(distmax))
prop.mat <-dist.daily_probs<-dist.propSingle<-matrix(nrow=locs,ncol=length(distmax))
for (j in 1:length(genvec)){
  for (n in 1:locs){
    # startSamp=n
    # startSamp<-sample(234,1,replace=T,prob=pop2019.town)
    startSamp<-sample(234,1,replace=T)
    
  for (k in 1:length(distmax)){
     ## across generations 
       prop.mat[n,k] <- sum(TranMatArray.1[startSamp,,genvec[j]][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
    ## single individual after infectious period   
       dist.propSingle[n,k]<-sum(tmpbase[startSamp,][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
  ## single individual daily 
       dist.daily_probs[n,k]<-sum(tmpbase_pre2[startSamp,][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
     ## random mobility  
       dist.prop.mat.rand[n,k] <-sum(pair.rand[startSamp,][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])
     ## raw facebook data
        dist.propFB[n,k]<-sum(cdr.mat.town[startSamp,][which(pairwise_geodist.town[startSamp,]>=distmin[k]&pairwise_geodist.town[startSamp,]<distmax[k])])

      }
  }
  # prop.mat.gens[j,]<-apply(prop.mat,2,function(x) weighted.mean(x,mrcaVec,na.rm=T ))
  prop.mat.gens[j,]<-apply(prop.mat,2,function(x) mean(x))
  
}

prop.mat.gens<-data.table(prop.mat.gens)
dist.prop.mat.list[[m]] <- prop.mat.gens
print(m)

###single person no generations infectious period
dist.propSingle.gens<-apply(dist.propSingle,2,function(x) mean(x))
dist.propSingle.gens<-data.table(t(dist.propSingle.gens))
dist.propSingle.list[[m]]<-dist.propSingle.gens


###single person no generations daily
dist.daily_probs.gens<-apply(dist.daily_probs,2,function(x) mean(x))
dist.daily_probs.gens<-data.table(t(dist.daily_probs.gens))
dist.daily_probs.list[[m]]<-dist.daily_probs.gens

}


### Find quantiles across bootstraps transmission event
dist.prop.mat2.tmp <- rbindlist(dist.prop.mat.list)
colnames(dist.prop.mat2.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
# dist.prop.mat2.tmp$gens <- rep(as.character(genvec),iters)
dist.prop.mat2.tmp$gens <- rep(as.character(genvec),nrow(pars))

### Find quantiles across bootstraps individual person across infectious period
dist.propSingle.tmp <- rbindlist(dist.propSingle.list)
colnames(dist.propSingle.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
# dist.prop.mat2.tmp$gens <- rep(as.character(genvec),iters)
dist.propSingle.tmp$gens <- "Individual_infec"

### Find quantiles across bootstraps individual person daily
dist.daily_probs.tmp <- rbindlist(dist.daily_probs.list)
colnames(dist.daily_probs.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
# dist.prop.mat2.tmp$gens <- rep(as.character(genvec),iters)
dist.daily_probs.tmp$gens <- "Individual_daily"

####random mobility
dist.prop.mat.rand.tmp<-data.table(dist.prop.mat.rand)
colnames(dist.prop.mat.rand.tmp) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
dist.prop.mat.rand.tmp$gens <- "Random"
####raw FB proportions
dist.propFB<-data.table(dist.propFB)
colnames(dist.propFB) <- c("0-10km","10-200km","200-500km","500-1000km",">1000km")
dist.propFB$gens <- "Meta"

dist.prop.mat2<-rbind(dist.prop.mat2.tmp,dist.prop.mat.rand.tmp,dist.propFB,dist.propSingle.tmp, dist.daily_probs.tmp)
dist.prop.mat2 <- melt(dist.prop.mat2)
colnames(dist.prop.mat2) <- c("gens","distances","value")
dist.prop.mat2$distances <- factor(dist.prop.mat2$distances,levels=c("0-10km","10-200km","200-500km","500-1000km",">1000km"),labels=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"))
genvec2<-c(1,10,"Individual_daily","Individual_infec","Random","Meta")
# genvec2<-c(1,50,"Random")



## Quantiles ### DOTS AND REDUCE THE NUMBER OF GENERATIONS
dists <- c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km")
dist.quants.list<-list()
for (i in 1:length(dists)) {
  dist.quants.tmp <- matrix(nrow=length(genvec2),ncol=6)
  for (j in 1:length(genvec2)){
    tmp <- subset(dist.prop.mat2, dist.prop.mat2$distances == dists[i] & dist.prop.mat2$gens==genvec2[j])
    dist.quants.tmp[j,1:3] <- as.numeric(quantile(tmp$value, probs=c(0.025,0.5,0.975)))
    dist.quants.tmp[j,4] <- as.numeric(mean(tmp$value))
    dist.quants.tmp[j,5] <- dists[i]
    dist.quants.tmp[j,6] <- genvec2[j]
  }
  dist.quants.tmp <- data.table(dist.quants.tmp)
  dist.quants.list[[i]] <- dist.quants.tmp
}
dist.quants <- rbindlist(dist.quants.list)
# colnames(dist.quants) <- c("median","lowerCI","upperCI","mean","distances","gens")
colnames(dist.quants) <- c("lowerCI","median","upperCI","mean","distances","gens")
dist.quants$lowerCI <- as.numeric(dist.quants$lowerCI)
dist.quants$median <- as.numeric(dist.quants$median)
dist.quants$mean <- as.numeric(dist.quants$mean)

dist.quants$upperCI <- as.numeric(dist.quants$upperCI)

# dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==10 | dist.quants$gens=="Random")
dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==10 | dist.quants$gens=="Random"
                           | dist.quants$gens=="Meta" | dist.quants$gens=="Individual_daily" | dist.quants$gens=="Individual_infec"  )

dist.quants.tmp3$lowerCI[dist.quants.tmp3$gens=="Random"] <- NA
dist.quants.tmp3$upperCI[dist.quants.tmp3$gens=="Random"] <- NA
dist.quants.tmp3$lowerCI[dist.quants.tmp3$gens=="Meta"] <- NA
dist.quants.tmp3$upperCI[dist.quants.tmp3$gens=="Meta"] <- NA



# dist.quants.tmp3.noPopWgt<-dist.quants.tmp3
# setwd("./ModelProjections/data/")
# save(dist.quants,file="dist.quants.posteriors.RData")
# save(dist.quants.tmp3.noPopWgt,file="dist.quants.posteriors.noPopwgt.RData")


##### save this file 17 May 2023
# save(dist.quants.tmp3,file="./ModelProjections/data/dist.quants.tmp3.posteriorsAll_adj.RData")

# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/proppergen/dist.quants.tmp3.posteriors.RData")

#####Combine two methods and plot
# 
ggplot(data=dist.quants.tmp3) +
  geom_point(data=dist.quants.tmp3,aes(x=distances,y=median,group=gens,shape=gens),size=3,position=position_dodge(width=0.35)) +
  geom_errorbar(data=dist.quants.tmp3,aes(ymin=lowerCI,ymax=upperCI,x=distances,group=gens),alpha=0.6,width=0.5,position=position_dodge(width=0.35))+
  xlab("Distance \n(from starting municipality)")+
  ylab("Proportion")+
  theme(panel.grid.minor  = element_blank(),)+
  theme_classic()+
  # scale_shape_manual(values=c(19,2,7),limits=c(1,10,"Random"), labels=c("1 Generation","10 Generations","Given random \nmobility")) +
  ylim(0,1)+
  theme(
    axis.text.x = element_text(angle=90,size=rel(2)), axis.text.y=element_text(size=rel(2)),
        # axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),legend.title = element_blank(),
        legend.position = c(.65,.75),legend.background = element_rect(color=NA, size=0.5, linetype="solid"))+

  # scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,50,"Random"), labels=c("1 Generation","50 Generations","Given random \nmobility")) +
  scale_x_discrete(limits=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"),labels=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"))+
  guides(fill="none")

# ggsave("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/propPerDist.10gens.parsUncertainty.pdf",width=4,height=6)
##################################



######################################################################################################
######PROVINCE LEVEL######
#######################################################################################################
# library(data.table)
# setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/RData_outputs/")
# load("pairwise_geodist.RData") # # [21_05_21_PariwiseDistance.R] ## Pairwise distance between towns in South Africa
# load("pop_2019.RData") 
# load("cdr.mat.one.RData")# # [Mobility_ManyMonthsSA.R] ## Probability of movement between each province normalized to carriage rates for each province
# cdr.mat<-cdr.mat.one
# 
# tn <- rownames(cdr.mat)
# pairwise_geodist<-pairwise_geodist[tn,tn]
# pop_2019 <-pop_2019[tn]
# pairwise_geodist[which(pairwise_geodist==0)]<-5
# nlocs=dim(pairwise_geodist)[1]
# distmin <- c(0,10,500,1000)
# distmax <- c(10,500,1000,1500)
# genvec <- c(1,10,50)
# # genvec <- c(1,6)
# 
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/TranMatArray.munic.RData")
# TranMatArray.1<-TranMatArray
# 
# TranMA.tmp <- array(data=NA, dim=c(nlocs,nlocs,length(genvec)))
# iters=1000
# pair.rand <- matrix(1/9,nrow=9,ncol=9)
# dist.prop.mat.list <-list()
# for (m in 1:iters){
#   dist.prop.mat <- matrix(nrow=length(genvec), ncol=length(distmax))
#   dist.prop.mat.rand <-vector(mode="numeric",length=length(distmax))
#   for (j in 1:length(genvec)){
#     prop.mat<-matrix(nrow=length(distmax),ncol=nlocs)
#     prop.mat.rand<-matrix(nrow=length(distmax),ncol=nlocs)
#     for (i in 1:nlocs){
#       # startSamp <- sample(234,1,replace=T,prob=mrcaVec)
#       startSamp <- sample(9,1,replace=T,prob=pop_2019)
#       
#       # startSamp <- sample(234,1,replace=T)
#       for (k in 1:length(distmax)){
#         prop.mat[k,i] <- sum(TranMatArray.1[startSamp,,genvec[j]][which(pairwise_geodist[startSamp,]>=distmin[k]&pairwise_geodist[startSamp,]<distmax[k])])
#         prop.mat.rand[k,i]<-sum(pair.rand[startSamp,][which(pairwise_geodist[startSamp,]>=distmin[k]&pairwise_geodist[startSamp,]<distmax[k])])
#       }
#       # print(i)
#     }
#     dist.prop.mat[j,] <- apply(prop.mat,1,function(x) mean(x))
#     dist.prop.mat.rand  <- apply(prop.mat.rand,1,function(x) mean(x))
#     # dist.prop.mat[j,] <- apply(prop.mat,1,function(x) weighted.mean(x,pop_2019))
#     # dist.prop.mat.rand  <- apply(prop.mat.rand,1,function(x) weighted.mean(x,pop_2019))
#   }
#   dist.prop.mat <- data.table(dist.prop.mat)
#   dist.prop.mat.list[[m]] <- dist.prop.mat
#   print(m)
# }
# 
# 
# 
# ### Find quantiles across bootstraps
# dist.prop.mat2.tmp <- rbindlist(dist.prop.mat.list)
# colnames(dist.prop.mat2.tmp) <- c("0-10km","10-500km","500-1000km","1000-1500km")
# dist.prop.mat2.tmp$gens <- rep(as.character(genvec),iters)
# dist.prop.mat.rand.tmp<-data.table(t(data.table(as.numeric(dist.prop.mat.rand))))
# colnames(dist.prop.mat.rand.tmp) <- c("0-10km","10-500km","500-1000km","1000-1500km")
# dist.prop.mat.rand.tmp$gens <- "Random"
# 
# dist.prop.mat2<-rbind(dist.prop.mat2.tmp,dist.prop.mat.rand.tmp)
# dist.prop.mat2 <- melt(dist.prop.mat2)
# colnames(dist.prop.mat2) <- c("gens","distances","value")
# dist.prop.mat2$distances <- factor(dist.prop.mat2$distances,levels=c("0-10km","10-500km","500-1000km","1000-1500km"),labels=c("Within Province","10-500km","500-1000km","1000-1500km"))
# genvec2<-c(1,10,50,"Random")
# # genvec2<-c(1,6,"Random")
# 
# ## Quantiles ### DOTS AND REDUCE THE NUMBER OF GENERATIONS
# dists <- c("Within Province","10-500km","500-1000km","1000-1500km")
# dist.quants.list<-list()
# for (i in 1:length(dists)) {
#   dist.quants.tmp <- matrix(nrow=length(genvec2),ncol=6)
#   for (j in 1:length(genvec2)){
#     tmp <- subset(dist.prop.mat2, dist.prop.mat2$distances == dists[i] & dist.prop.mat2$gens==genvec2[j])
#     dist.quants.tmp[j,1:3] <- as.numeric(quantile(tmp$value, probs=c(0.025,0.5,0.975)))
#     dist.quants.tmp[j,4] <- as.numeric(mean(tmp$value))
#     dist.quants.tmp[j,5] <- dists[i]
#     dist.quants.tmp[j,6] <- genvec2[j]
#   }
#   dist.quants.tmp <- data.table(dist.quants.tmp)
#   dist.quants.list[[i]] <- dist.quants.tmp
# }
# dist.quants <- rbindlist(dist.quants.list)
# colnames(dist.quants) <- c("lowerCI","median","upperCI","mean","distances","gens")
# dist.quants$lowerCI <- as.numeric(dist.quants$lowerCI)
# dist.quants$median <- as.numeric(dist.quants$median)
# dist.quants$mean <- as.numeric(dist.quants$mean)
# 
# dist.quants$upperCI <- as.numeric(dist.quants$upperCI)
# 
# # dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==6 | dist.quants$gens=="Random")
# dist.quants.tmp3 <- subset(dist.quants, dist.quants$gens==1 | dist.quants$gens==50 | dist.quants$gens=="Random")
# 
# dist.quants.tmp3$lowerCI[dist.quants.tmp3$gens=="Random"] <- NA
# dist.quants.tmp3$upperCI[dist.quants.tmp3$gens=="Random"] <- NA
# 
# ggplot(data=dist.quants.tmp3) +
#   geom_point(data=dist.quants.tmp3,aes(x=distances,y=median,shape=gens),size=3,position=position_dodge(width=0.35)) +
#   geom_errorbar(aes(ymin=lowerCI,ymax=upperCI,shape=gens,x=distances),alpha=0.6,width=0.2,position=position_dodge(width=0.35))+
#   xlab("Distance \n(from starting province)")+
#   ylab("Proportion")+
#   theme(panel.grid.minor  = element_blank(),)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90,size=18), axis.text.y=element_text(size=18),
#         axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
#         legend.text = element_text(size=12),legend.title = element_text(size=12),
#         legend.position = c(.8,.8),legend.background = element_rect(color="black", size=0.5, linetype="solid"))+
#   scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,50,"Random"), labels=c("1 Generation","2 years","Given random \nmobility")) +
#   # scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,10,"Random"), labels=c("1 Generation","10 Generations","Given random \nmobility")) +
#   scale_x_discrete(limits=c("Within Province","10-500km","500-1000km","1000-1500km"),labels=c("Within Province","10-500km","500-1000km","1000-1500km"))+
#   guides(fill="none") +
#   ylim(0,1)
# 

