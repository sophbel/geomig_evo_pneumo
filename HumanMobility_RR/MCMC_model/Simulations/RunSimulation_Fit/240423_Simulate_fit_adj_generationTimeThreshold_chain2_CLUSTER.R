setwd("/data/pam/team284/sb62/scratch/Migration/SouthAfrica/mobility_model/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
####Fit model with simulation
iters=20000
scale=.07
boot=2
#########Set up simulation
library(dplyr)
library(data.table)
library(ggplot2)
library(abind)
library(doParallel)
library(ucminf)
library(doMC)
library(Rcpp)
library(RcppEigen)
library(Rfast)
library(coda)
library(fmcmc)
load('./modelinput_data/cdr.mat.one.RData')
cdr.mat<-cdr.mat.one
load('./modelinput_data/cdr.mat.town.one.RData')
cdr.mat.town<-cdr.mat.town.one
load('./modelinput_data/pairwise_geodist.RData')
load("./modelinput_data/dat.tmp.allser.RData")
load('./modelinput_data/pop_2019.RData')
load("./modelinput_data/pop_municipality.2017LS.RData")
sourceCpp("./MCMC_model/MatrixMultiplication.cpp")

#####index for province and municipality match
nloc.munic<-length(pop2019.town)
splitnames <- matrix(nrow=nloc.munic,ncol=2)
for (us in 1:nloc.munic) {
  for (prov in 1:2) { 
    splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
  }
}
Eigencpp=TRUE
a<-match(colnames(cdr.mat),names(pop_2019))
pop_2019<-pop_2019[a]

a<-match(colnames(cdr.mat.town),names(pop2019.town))
pop2019.town<-pop2019.town[a]


## Simulation parameters
nsims<-1  ## Number of sepaarte simulations to run
burnin<-20
# nGen2<-nGen<-1:30
# pGen2<-pGen<-rep(1/30,30)
nGen2<-nGen<-1:120
pGen2<-pGen<-rep(1/120,120)
ntimes<-50000 ## No. of pairs per simulation per observation

### True sampling from true data set.
dat.inMaster<-do.call("rbind",dat.tmp.allser)
probByloc<-hist(c(dat.inMaster[,3],dat.inMaster[,4]),breaks=1:10-0.5,plot=F)$counts
probSample<-probByloc
# probSample<-probSample/max(probSample)
probSample<-probSample/sum(probSample)

##### Add parameter to Mobility data 
parHome<--2
min.range=-0.04
max.range=0.6
nlocs=nloc=9
npop<-pop_2019
pHome<-npop/sum(npop)


#########MUNICIPALITY LEVEL
# tmp.pHome.PROV<-extTranMatDat.tmp$popbyCell
# tmp.pHome.PROV<-tmp.pHome.PROV/sum(tmp.pHome.PROV)

### Create mobility matrix for susceptible individuals
tmpbase_pre<-cdr.mat.town
# tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
tmppar1 <- exp(parHome)/(1+exp(parHome))
tmppar <- min.range+tmppar1*(max.range-min.range)
tmpdiag<-diag(tmpbase_pre)-tmppar
tmpdiag[which(tmpdiag>0.99999)]<-0.99999
diag(tmpbase_pre)<-0
tmpbase_pre1<-sweep(tmpbase_pre,1,rowSums(tmpbase_pre),"/")
tmpbase_pre2<-sweep(tmpbase_pre1,1,(1-tmpdiag)/(1-diag(tmpbase_pre1)),"*")
diag(tmpbase_pre2)<-tmpdiag



#### adjust so it is mobility across the infectious period
timeWindow<-35
probStay<-1-(diag(tmpbase_pre2))^timeWindow
tmp<-tmpbase_pre2
diag(tmp)<-0
tmp1<-sweep(tmp,1,rowSums(tmp),"/")
tmp2<-sweep(tmp1,1,(1-probStay)/(1-diag(tmp1)),"*")
diag(tmp2)<-probStay
tmpbase<-tmp2

####### Collapse mobility matrix to 9X9 here
nloc.munic=nrow(tmpbase)

splitnames <- matrix(nrow=nloc.munic,ncol=2)
for (us in 1:nloc.munic) {
  for (prov in 1:2) { 
    splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
  }
}

tmpbase.1 <- apply(tmpbase, 1, function(x){
  split(x,splitnames[,2]) %>% sapply(sum)
})  %>% transpose()

tmpbase2 <- apply(tmpbase.1,2,function(x) {
  mapply( weighted.mean
          , x = split(x,splitnames[,2])
          , w = split(pop2019.town,splitnames[,2])
  )
})

move3<-tcrossprod(tmpbase2,tmpbase2)
move4<-sweep(move3,2,pHome,"*")
TranMat<-sweep(move4,1,rowSums(move4),"/")
##############

#SIMULATION FUNCTION
SimulationFunction<-function(){
  
  start=sample(nloc,ntimes,prob=pHome,replace=T)
  for (jj in 1:burnin){
    for (kk in 1:ntimes){
      start[kk]=sample(nloc,1,prob=TranMat[start[kk],])
    }
  }
  whereEnd2<-whereEnd<-start
  nogens<-sample(nGen,ntimes,replace=T,prob=pGen)
  nogens2<-sample(nGen2,ntimes,replace=T,prob=pGen2)
  for (kk in 1:ntimes){
    for (jj in 1:nogens[kk]){
      whereEnd[kk]=sample(nloc,1,prob=TranMat[whereEnd[kk],])
    }
    for (jj in 1:nogens2[kk]){
      whereEnd2[kk]=sample(nloc,1,prob=TranMat[whereEnd2[kk],])
    }
  }
  
  output<-list(start, whereEnd,whereEnd2,nogens,nogens2)
  return(output)
}

# for(boot in 1:3){
print("Running Simulation:")
  print(boot)
sim<-SimulationFunction()

start<-sim[[1]]
whereEnd<-sim[[2]]
whereEnd2<-sim[[3]]
nogens<-sim[[4]]
nogens2<-sim[[5]]


##Prepare data
dat.in.all<-as.data.frame(cbind(id1=1:length(start),id2=1:length(start),
                                start=start,loc1=whereEnd2,loc2=whereEnd,
                                nogens=nogens,nogens2=nogens2,meanGen=nogens+nogens2))

a<-as.matrix(dat.in.all[,4:5])
dat.in.all$dist<-pairwise_geodist[a]
a<-as.matrix(dat.in.all[,c(3,5)])
dat.in.all$distStart1<-pairwise_geodist[a]
a<-as.matrix(dat.in.all[,c(3,4)])
dat.in.all$distStart2<-pairwise_geodist[a]



####SUB-Sample those in the same place
# samp1<-rbinom(nrow(dat.in.all),1,prob=probSample[dat.in.all[,"loc1"]])
# samp2<-rbinom(nrow(dat.in.all),1,prob=probSample[dat.in.all[,"loc2"]])
# dat.in.tmp<-dat.in.all[which(samp1==1&samp2==1),]


### Randomly sub-sample
# randsamp<-rep(1/9,9)
# samp1<-rbinom(nrow(dat.in.all),1,prob=randsamp[dat.in.all[,"loc1"]])
# samp2<-rbinom(nrow(dat.in.all),1,prob=randsamp[dat.in.all[,"loc2"]])
# dat.in.tmp<-dat.in.all[which(samp1==1&samp2==1),]
####or
samp<-sample(nrow(dat.in.all),5000,replace=F)
dat.in.tmp<-dat.in.all[samp,]


# ###Sample 3X more 
# samp1<-rbinom(nrow(dat.in.all)*5,1,prob=probSample[dat.in.all[,"loc1"]])
# samp2<-rbinom(nrow(dat.in.all)*5,1,prob=probSample[dat.in.all[,"loc2"]])
# samp.tmp<-which(samp1==1&samp2==1)
# samp.tmp<-samp.tmp[which(samp.tmp<50000)]
# dat.in.tmp<-dat.in.all[samp.tmp,]

print("Dimensions sub-sample")
print(dim(dat.in.tmp))
locs1<-table(dat.in.tmp$loc1)
locs2<-table(dat.in.tmp$loc2)
print("Number per province")
print(locs1)
print(locs2)


###################
##############REFIT MODEL USING SIMULATION##########################################
calcAllProbs=FALSE
singPar=TRUE

ncore=1
# min.range<-(-0.04)
# max.range<-0.6

extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999

source("./MCMC_model/Simulations/Likelihoods/LikFunc.mcmc.SIMS.Munic_adj.R")
# source("./MCMC_model/Simulations/Likelihoods/LikFunc.mcmc.SIMS.Munic_adj_singlePar.R")

nInfecLoc<-table(dat.in.all$start)

if(singPar==TRUE){
  startPar<-runif(1,-3,0)
  
}else{
  startPar<-c(-2.5,rep(0,8))
  
}

#############################################
##########set the generation times to 15
max.no.gens<-15 ### 
maxGen<-max.no.gens#### 
dat.in2<-subset(dat.in.tmp,dat.in.tmp$nogens<=max.no.gens&dat.in.tmp$nogens2<=max.no.gens)

##Probability of being in each generation set at true generations
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}

npairs=nrow(dat.in2)
print(paste0("Number of pairs in subsample=",npairs))

#######RUN MCMC#####
print("Running MCMC")
# startPar<-c(-2.5,rep(0,8))
start.time<-Sys.time()
sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = scale),progress = interactive())
end.time<-Sys.time()
print(end.time-start.time)
if(singPar==TRUE){
  save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres15_singPar.RData"))
  save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres15_singPar.RData"))
  save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres15_singPar.RData"))###SAVE OUTPUT
}else{
save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres15.RData"))
save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres15.RData"))
save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres15.RData"))###SAVE OUTPUT
}
######################################################
#########set the generation times to 30
max.no.gens<-30 ### 
maxGen<-max.no.gens#### 
dat.in2<-subset(dat.in.tmp,dat.in.tmp$nogens<=max.no.gens&dat.in.tmp$nogens2<=max.no.gens)

####Probability of being in each generation set at true generations
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}

npairs=nrow(dat.in2)
print(paste0("Number of pairs in subsample=",npairs))

#######RUN MCMC#######################
print("Running MCMC")
# startPar<-c(-2.5,rep(0,8))
start.time<-Sys.time()
sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = scale),progress = interactive())
end.time<-Sys.time()
print(end.time-start.time)
if(singPar==TRUE){
  save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres30_singPar.RData"))
  save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres30_singPar.RData"))
  save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres30_singPar.RData"))###SAVE OUTPUT
}else{
save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres30.RData"))
save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres30.RData"))
save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres30.RData"))###SAVE OUTPUT
}
########################################################################
#########set the generation times to 60
max.no.gens<-60 ### 
maxGen<-max.no.gens#### 
dat.in2<-subset(dat.in.tmp,dat.in.tmp$nogens<=max.no.gens&dat.in.tmp$nogens2<=max.no.gens)

####Probability of being in each generation set at true generations
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}

npairs=nrow(dat.in2)
print(paste0("Number of pairs in subsample=",npairs))

#######RUN MCMC#######################
print("Running MCMC")
# startPar<-c(-2.5,rep(0,8))
# startPar<-runif(1,-3,0)

start.time<-Sys.time()
sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = scale),progress = interactive())
end.time<-Sys.time()
print(end.time-start.time)
if(singPar==TRUE){
  save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres60_singPar.RData"))
  save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres60_singPar.RData"))
  save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres60_singPar.RData"))###SAVE OUTPUT
}else{
save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres60.RData"))
save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres60.RData"))
save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres60.RData"))###SAVE OUTPUT
}

#########set the generation times to 90
max.no.gens<-90 ### 
maxGen<-max.no.gens#### 
dat.in2<-subset(dat.in.tmp,dat.in.tmp$nogens<=max.no.gens&dat.in.tmp$nogens2<=max.no.gens)

####Probability of being in each generation set at true generations
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}

npairs=nrow(dat.in2)
print(paste0("Number of pairs in subsample=",npairs))

#######RUN MCMC#######################
print("Running MCMC")
# startPar<-c(-2.5,rep(0,8))
# startPar<-runif(1,-3,0)

start.time<-Sys.time()
sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = scale),progress = interactive())
end.time<-Sys.time()
print(end.time-start.time)
if(singPar==TRUE){
  save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres90_singPar.RData"))
  save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres90_singPar.RData"))
  save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres90_singPar.RData"))###SAVE OUTPUT
}else{
save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres90.RData"))
save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres90.RData"))
save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres90.RData"))###SAVE OUTPUT
}

#########set the generation times to 120
max.no.gens<-120 ### 
maxGen<-max.no.gens#### 
dat.in2<-subset(dat.in.tmp,dat.in.tmp$nogens<=max.no.gens&dat.in.tmp$nogens2<=max.no.gens)

####Probability of being in each generation set at true generations
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}

npairs=nrow(dat.in2)
print(paste0("Number of pairs in subsample=",npairs))

#######RUN MCMC#######################
print("Running MCMC")


start.time<-Sys.time()
sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = scale),progress = interactive())
end.time<-Sys.time()
print(end.time-start.time)
if(singPar==TRUE){
  save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres120_singPar.RData"))
  save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres120_singPar.RData"))
  save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres120_singPar.RData"))###SAVE OUTPUT
}else{
save(dat.in2,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in2",scale,".",iters,".",boot,"_adj_genthres120.RData"))
save(dat.in.all,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/dat.in.all",scale,".",iters,".",boot,"_adj_genthres120.RData"))
save(sim.ans,file=paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.",scale,".",iters,".",boot,"_adj_genthres120.RData"))###SAVE OUTPUT
}
# }
###############################################################################
