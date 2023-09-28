setwd( "/data/pam/team284/sb62/scratch/Migration/SouthAfrica/mobility_model/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
### Municipality Level Function
iters=20000
# for (chain in 1:3){
chain=1
cluster=FALSE
##### Adding variables from other scripts and loading libraries
library(data.table)
library(doParallel)
library(ucminf)
library(doMC)
library(Rcpp)
library(RcppEigen)
library(Rfast)
library(abind)
library(ggplot2)
library(data.table)
library(fmcmc)
library(coda)
library(dplyr)
#### Bring in data
load('./modelinput_data/cdr.mat.one.RData')
cdr.mat<-cdr.mat.one
# load('./modelinput_data/cdr.mat.IP.RData')
# cdr.mat<-cdr.mat.IP
load('./modelinput_data/pop_2019.RData')
load('./modelinput_data/dat.tmp.allser.RData')
# load('./modelinput_data/tMRCAs.RData')
load('./modelinput_data/pairwise_geodist.RData')
load("./modelinput_data/pop_municipality.2017LS.RData") 
load("./modelinput_data/pairwise_geodist.town.RData") 
load('./modelinput_data/cdr.mat.town.one.RData')
cdr.mat.town<-cdr.mat.town.one
# load('./modelinput_data/cdr.mat.munic.IP.RData')
# cdr.mat.town<-cdr.mat.munic.IP
###
nloc=nlocs=9
nloc.munic=234
maxGen<-100

maxTranGens=maxGen

## Gamma distribution for a generation time of 35 days- same as lines above
genTime<- 35
varGen<-(genTime)^2
shape=genTime^2/varGen
scale=varGen/genTime

# max.no.gens<-round(3650/genTime,0)
max.no.gens<-maxGen

plot(1:1000,dgamma(1:1000,shape=shape,scale=scale))



pFunc<-function(x){
  a<-dgamma(x*365,shape=shape*1:1000,scale=scale)
  den<-sum(a,na.rm=T)
  b<-a[1:max.no.gens]/den
  return(b)
}


prob.tMRCA.out<-list()
for (ser in 1:9){
  gD<-tMRCAs[[ser]]
  prob.gdist<-array(0,c(nrow(gD),nrow(gD),max.no.gens))
  for (i in 1:nrow(gD)){
    for (j in 1:nrow(gD)){
      if(i==j)next
      if((gD[i,j]*365)>5000)next#### This results in an error if 5000 is too small or too large
      prob.gdist[i,j,]<-pFunc(x=gD[i,j])
    }
  }
  attr(prob.gdist,"data")<-attr(gD,"data")
  prob.tMRCA.out[[ser]]<-prob.gdist
  print(ser)
}

probByGen.tmp=prob.tMRCA.out
save(prob.tMRCA.out,file="./modelinput_data/prob.tMRCA.out.RData")

## set provs order
provs<-colnames(cdr.mat)

### Adding variables
Eigencpp=TRUE
calcAllProbs=FALSE

maxGen=maxTranGens

sourceCpp("./MCMC_model/MatrixMultiplication.cpp")

## original mobility
# source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic.R")

## adjusting the mobility to span the infectious period rather than being daily
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_adj.R")


nam.a<-colnames(cdr.mat)
pop_2019<-pop_2019[nam.a]
pairwise_geodist<-pairwise_geodist[nam.a,nam.a]

nam.b<-colnames(cdr.mat.town)
pop2019.town<-pop2019.town[nam.b]
pairwise_geodist.town<-pairwise_geodist.town[nam.b,nam.b]
### Prepare data
nByGPSC<-rep(NaN,9)
for (i in 1:9)nByGPSC[i]<-nrow(dat.tmp.allser[[i]])
dat.inMaster<-do.call("rbind",dat.tmp.allser)
dat.inMaster<-dat.inMaster[which(is.na(dat.inMaster$totTimeDays)==F),]
dat.inMaster$GPSC<-rep(1:9,nByGPSC)

dat.inMaster$year1<-floor(dat.inMaster$cy1)
dat.inMaster$year2<-floor(dat.inMaster$cy2)
a<-as.matrix(dat.inMaster[,3:4])
dat.inMaster$SDist<-pairwise_geodist[a]


### Sampling probability by year, GPSC and province
NoDetectByYearGPSC<-table(c(dat.inMaster$loc1,dat.inMaster$loc2),c(dat.inMaster$GPSC,dat.inMaster$GPSC),c(dat.inMaster$year1,dat.inMaster$year2))
un.years<-sort(unique(dat.inMaster$year1))

avNoDetectByloc<-table(c(dat.inMaster$loc1,dat.inMaster$loc2))
dat.inMaster$YearREF1<-match(dat.inMaster$year1,un.years)
dat.inMaster$YearREF2<-match(dat.inMaster$year2,un.years)

dat.inMaster$SDist[which(dat.inMaster$SDist==0)]<-10
pairwise_geodist[which(pairwise_geodist==0)]<-10



### Analysis dataset
ncore=1
min.range<-(-0.04)
max.range<-0.6

extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999

dat.in2<-dat.inMaster[which(dat.inMaster$totTimeDays<3500),]
npairs=nrow(dat.in2)

### Run MCMC
par1<-runif(1,-3,0)
par2_8<-runif(8,0,0.9999)
startPar<-c(par1,par2_8)
ans.munic <- MCMC(likFunc.munic,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = .08),thin=5)
save(ans.munic,file=paste0("./MCMC_model/outputs/mechan_model/ans.munic",chain,".",iters,".08_adj",".RData"))





