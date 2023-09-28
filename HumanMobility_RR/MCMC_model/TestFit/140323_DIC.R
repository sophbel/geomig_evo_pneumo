setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
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
###### Data for recalculating the likelihood########
#### Bring in data
load('./modelinput_data/cdr.mat.one.RData')
cdr.mat<-cdr.mat.one
load('./modelinput_data/pop_2019.RData')
load('./modelinput_data/dat.tmp.allser.RData')
load('./modelinput_data/tMRCAs.RData')
load('./modelinput_data/pairwise_geodist.RData')
load("./modelinput_data/pop_municipality.2017LS.RData") 
load("./modelinput_data/pairwise_geodist.town.RData") 
load('./modelinput_data/cdr.mat.town.one.RData')
cdr.mat.town<-cdr.mat.town.one

########### Set up data
nloc=nlocs=9
nloc.munic=234
maxGen<-100

maxTranGens=maxGen

## Gamma distribution for a generation time of 35 days- same as lines above
genTime<- 35
varGen<-(genTime)^2
shape=genTime^2/varGen
scale=varGen/genTime

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
## set provs order
provs<-colnames(cdr.mat)

### Adding variables
Eigencpp=TRUE
calcAllProbs=FALSE

maxGen=maxTranGens

sourceCpp("./MCMC_model/MatrixMultiplication.cpp")

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

##### for analysis
ncore=1
min.range<-(-0.04)
max.range<-0.6

extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999

dat.in2<-dat.inMaster[which(dat.inMaster$totTimeDays<3500),]
npairs=nrow(dat.in2)

##########Calculate the DIC########
### mean of each parameter
par_means<-function(chains){
  mean_pars<-colMeans(chains)
  return(mean_pars)
}

#--- DIC Calculation 
library(EpiILM)


calc_dic <- function(data, chains){


  # log-likelihood at each iteration (Dbar)
  liki <- data
  # print(liki)
  # log-likelihood at mean of posterior
  mpars<-par_means(chains)
  plik <-  likFunc.munic(mpars)

  ## mean likelihood 
  menliki<-mean(liki)

  ## MRC Version 2002
  D.bar = -2 * mean(liki)
  D.theta_bar = -2 * plik
  pD = D.bar - D.theta_bar
  DIC = pD + D.bar
  dic_vec<-c(DIC,pD,D.bar)
  names(dic_vec)<-c("dic","pD","DBar")
  # 
  ## Gelman Version
  # D.bar = -2 * mean(liki)
  # coliki<-melt(liki)$value
  # D.theta = -2 * coliki
  # pD = var(D.theta) / 2
  # DIC = pD + D.bar
  # dic_vec<-c(DIC,pD,D.bar)
  # names(dic_vec)<-c("dic","pD","DBar")
  
  ## Megan old vers
  # Pd <- 2*menliki - 2*plik
  # dic <- -2*plik + 2*Pd
  # dic_vec<-c(dic,menliki,plik,Pd)
  # names(dic_vec)<-c("dic","lik_chain_mean","plik","Pd")
  
  return(dic_vec)
}
 
# library(loo)
# waic(log_liks)
# 
# ##### Set wd#####
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
dic_mat<-matrix(nrow=5,ncol=4)
######################################################################
############## Human Mobility Model####
##### Model outputs ######## ##########################################
load("./MCMC_model/outputs/mechan_model/ans.munic1.20000.08_adj.RData")
chain1<-ans.munic
chain1<-chain1[3000:4000,]
load("./MCMC_model/outputs/mechan_model/ans.munic2.20000.08_adj.RData")
chain2<-ans.munic
chain2<-chain2[3000:4000,]
load("./MCMC_model/outputs/mechan_model/ans.munic3.20000.08_adj.RData")
chain3<-ans.munic
chain3<-chain3[3000:4000,]
chains<-rbind(chain1,chain2,chain3)

#### LLs from runs #####
LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain1_adjhum.o_LLs.txt",header=FALSE)
LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs1<-LLs1[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain2_adjhum.o_LLs.txt",header=FALSE)
LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs2<-LLs2[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain3_adjhum.o_LLs.txt",header=FALSE)
LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs3<-LLs3[3000:4000]

##### bind all likelihoods from each iteration in each chain
log_liks<-cbind(LLs1,LLs2,LLs3)
#--- likelihood function to apply across all
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_adj.R")

mpars<-par_means(chains)
plik <-  likFunc.munic(mpars)
# epidic(3000, 4000, log_liks, plik)
dic<-calc_dic(log_liks,chains)
dic<-c(dic,"hum_model")
dic_mat[1,]<-dic
# 

#
# ######################################################################
# ############## Gravity Model JUST GAMMA####
# ##### Model outputs ######## ##########################################
### all chains
iters=20000
load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic1.",iters,".08_gravity_adj_g.RData"))
chain1<-ans.munic
chain1<-chain1[3000:4000,]
load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic2.",iters,".08_gravity_adj_g.RData"))
chain2<-ans.munic
chain2<-chain2[3000:4000,]
load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic4.",iters,".08_gravity_adj_g.RData"))
chain3<-ans.munic
chain3<-chain3[3000:4000,]

chains<-rbind(chain1,chain2,chain3)
# chains<-rbind(chain1,chain3)



#### LLs from runs #####
LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain7_g_g.o_LLs.txt",header=FALSE)
LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs1<-LLs1[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain8_g_g.o_LLs.txt",header=FALSE)
LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs2<-LLs2[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain10_g_g.o_LLs.txt",header=FALSE)
LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs3<-LLs3[3000:4000]


##### bind all likelihoods from each iteration in each chain
log_liks<-cbind(LLs1,LLs2,LLs3)
# log_liks<-cbind(LLs1,LLs3)

#--- likelihood function to apply across all
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_g.R")
dic2<-calc_dic(log_liks,chains)
dic2<-c(dic2,"g_model")
dic_mat[2,]<-dic2


# ######################################################################
# ############## Gravity Model BETA GAMMA####
# ##### Model outputs ######## ##########################################
##### Model outputs ######## ##########################################
# load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic1.20000.05_gravity_adj_bg_scale.RData")
# chain1<-ans.munic
# chain1<-chain1[3000:4000,]
load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic2.20000.05_gravity_adj_bg_scale.RData")
chain2<-ans.munic
chain2<-chain2[3000:4000,]
load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic3.20000.06_gravity_adj_bg_scale.RData")
chain3<-ans.munic
chain3<-chain3[3000:4000,]

# chains<-rbind(chain1,chain2,chain3)
chains<-rbind(chain2,chain3)




#### LLs from runs #####
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/runs2/chain1_bg06.o_LLs.txt",header=FALSE)
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain21_bg20k.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs1<-LLs1[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/runs2/chain2_bg052.o_LLs.txt",header=FALSE)
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain22_bg20k.o_LLs.txt",header=FALSE)
LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs2<-LLs2[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/runs2/chain3_bg052.o_LLs.txt",header=FALSE)
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain23_bg20k.o_LLs.txt",header=FALSE)
LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs3<-LLs3[3000:4000]


### bind all likelihoods from each iteration in each chain
# log_liks<-cbind(LLs1,LLs2,LLs3)
log_liks<-cbind(LLs2,LLs3)

# log_liks<-LLs1

#--- likelihood function to apply across all
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_bg.R")
dic3<-calc_dic(log_liks,chains)
dic3<-c(dic3,"bg_model")

dic_mat[3,]<-dic3

####overall
dic_mat<-data.table(dic_mat)
colnames(dic_mat)<-c("DIC","pD","D_bar","model")
dic_mat<-data.frame(dic_mat)
col.nums<-c(1:3)
dic_mat[col.nums]<-sapply(dic_mat[col.nums],as.numeric)
dic_mat<-dic_mat[1:3,]
dic_mat$DIC_diff<-dic_mat$DIC-min(dic_mat$DIC)


ggplot(dic_mat)+
  geom_point(aes(x=model,y=DIC))


# ######################################################################
# ############## Gravity Model GAMMA + Parameter####
# ##### Model outputs ######## ##########################################
load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic1.",iters,".08_gravity_adj_gPar_meta2.RData"))
chain1<-ans.munic
load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic2.",iters,".08_gravity_adj_gPar_meta2.RData"))
chain2<-ans.munic
load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic3.",iters,".08_gravity_adj_gPar_meta2.RData"))
chain3<-ans.munic
chains<-rbind(chain1,chain2,chain3)
# chains<-rbind(chain1,chain3)
#### LLs from runs #####
LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain10_g_gPar.o_LLs.txt",header=FALSE)
LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs1<-LLs1[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain11_g_gPar.o_LLs.txt",header=FALSE)
LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs2<-LLs2[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain12_g_gPar.o_LLs.txt",header=FALSE)
LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs3<-LLs3[3000:4000]
##### bind all likelihoods from each iteration in each chain
log_liks<-cbind(LLs1,LLs2,LLs3)
# log_liks<-cbind(LLs1,LLs3)
# log_liks<-LLs1

#--- likelihood function to apply across all
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_gPar.R")
dic4<-calc_dic(log_liks,chains)
dic4<-c(dic4,"gPar_model")

dic_mat[4,]<-dic4


# ######################################################################
# ############## Gravity Model BETA GAMMA + Parameter####
# ##### Model outputs ######## ##########################################
# load(paste0("./MCMC_model/outputs/gravity_model/grav_bgPar/ans.munic1.",iters,".08_gravity_adj_gPar_meta2.RData"))
# chain1<-ans.munic
load(paste0("./MCMC_model/outputs/gravity_model/grav_bgPar/ans.munic2.",iters,".08_gravity_adj_bgPar_meta.RData"))
chain2<-ans.munic
load(paste0("./MCMC_model/outputs/gravity_model/grav_bgPar/ans.munic3.",iters,".08_gravity_adj_bgPar_meta.RData"))
chain3<-ans.munic
# chains<-rbind(chain1,chain2,chain3)
chains<-rbind(chain2,chain3)
#### LLs from runs #####
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bgPar/runs/chain1_g_bgPar.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs1<-LLs1[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_bgPar/runs/chain2_bgPar.o_LLs.txt",header=FALSE)
LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs2<-LLs2[3000:4000]
LLs_all<-read.csv("./MCMC_model/RunModel/grav_bgPar/runs/chain3_bgPar.o_LLs.txt",header=FALSE)
LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
LLs3<-LLs3[3000:4000]
##### bind all likelihoods from each iteration in each chain
# log_liks<-cbind(LLs1,LLs2,LLs3)
log_liks<-cbind(LLs2,LLs3)
# log_liks<-LLs1

#--- likelihood function to apply across all
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_bgPar.R")
dic5<-calc_dic(log_liks,chains)
dic5<-c(dic5,"bgPar_model")

dic_mat[5,]<-dic5













# 
# 
# 
# 
# 
# 
# 
# ###################
# ###Test with EpiLMC package#######
# par_means<-function(chains){
#   mean_pars<-colMeans(chains)
#   return(mean_pars)
# }
# 
# #--- DIC Calculation 
# library(EpiILM)
# setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
# dic<-vector(mode="numeric",4)
# ############## Human Mobility Model####
# load("./MCMC_model/outputs/mechan_model/ans.munic1.20000.08_adj.RData")
# chain1<-ans.munic
# load("./MCMC_model/outputs/mechan_model/ans.munic2.20000.08_adj.RData")
# chain2<-ans.munic
# load("./MCMC_model/outputs/mechan_model/ans.munic3.20000.08_adj.RData")
# chain3<-ans.munic
# chains<-rbind(chain1,chain2,chain3)
# #### LLs from runs
# LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain1_adjhum.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain2_adjhum.o_LLs.txt",header=FALSE)
# LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/HumMob/runs/chain3_adjhum.o_LLs.txt",header=FALSE)
# LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# ##### bind all likelihoods from each iteration in each chain
# log_liks<-cbind(LLs1,LLs2,LLs3)
# #--- likelihood function to apply across all
# source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_adj.R")
# mpars<-par_means(chains)
# plik <-  likFunc.munic(mpars)
# dic[1]<-epidic(3000, 4000, log_liks, plik)
# 
# ####### Gravity model bg############
# load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic1.20000.08_gravity_adj_bg.RData")
# chain1<-ans.munic
# load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic2.20000.08_gravity_adj_bg.RData")
# chain2<-ans.munic
# load("./MCMC_model/outputs/gravity_model/grav_bg/ans.munic3.20000.08_gravity_adj_bg.RData")
# chain3<-ans.munic
# chains<-rbind(chain1,chain2,chain3)
# #### LLs from runs
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain21_bg20k.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain22_bg20k.o_LLs.txt",header=FALSE)
# LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_bg/runs/chain23_bg20k.o_LLs.txt",header=FALSE)
# LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# log_liks<-LLs1
# source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_bg.R")
# mpars<-par_means(chains)
# plik <-  likFunc.munic(mpars)
# dic[2]<-epidic(3000, 4000, log_liks, plik)
# 
# ####### Gravity model g############
# iters=20000
# load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic1.",iters,".08_gravity_adj_g.RData"))
# chain1<-ans.munic
# load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic2.",iters,".08_gravity_adj_g.RData"))
# chain2<-ans.munic
# load(paste0("./MCMC_model/outputs/gravity_model/grav_g/ans.munic4.",iters,".08_gravity_adj_g.RData"))
# chain3<-ans.munic
# chains<-rbind(chain1,chain2,chain3)
# #### LLs from runs 
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain7_g_g.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain8_g_g.o_LLs.txt",header=FALSE)
# LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_g/runs/chain10_g_g.o_LLs.txt",header=FALSE)
# LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# ##### bind all likelihoods from each iteration in each chain
# log_liks<-cbind(LLs1,LLs2,LLs3)
# #--- likelihood function to apply across all
# source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_g.R")
# mpars<-par_means(chains)
# plik <-  likFunc.munic(mpars)
# dic[3]<-epidic(3000, 4000, log_liks, plik)
# 
# 
# ####### Gravity model g with meta###########
# iters=20000
# load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic1.",iters,".08_gravity_adj_gPar_meta2.RData"))
# chain1<-ans.munic
# # load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic2.",iters,".08_gravity_adj_gPar_meta2.RData"))
# # chain2<-ans.munic
# load(paste0("./MCMC_model/outputs/gravity_model/grav_gPar/ans.munic3.",iters,".08_gravity_adj_gPar_meta2.RData"))
# chain3<-ans.munic
# chains<-rbind(chain1,chain2,chain3)
# #### LLs from runs 
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain10_g_gPar.o_LLs.txt",header=FALSE)
# LLs1 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain11_g_gPar.o_LLs.txt",header=FALSE)
# LLs2 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# LLs_all<-read.csv("./MCMC_model/RunModel/grav_gPar/runs/chain12_g_gPar.o_LLs.txt",header=FALSE)
# LLs3 <- LLs_all$V1[seq(1, length(LLs_all$V1), 5)]
# ##### bind all likelihoods from each iteration in each chain
# # log_liks<-cbind(LLs1,LLs2,LLs3)
# log_liks<-cbind(LLs1,LLs3)
# 
# #--- likelihood function to apply across all
# source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_grav_gPar.R")
# mpars<-par_means(chains)
# plik <-  likFunc.munic(mpars)
# lls<-melt(log_liks)$value
# dic[4]<-epidic(3000, 4000, lls, plik)
# 
# names(dic)<-c("human_mob","grav_bg","grav_g","grav_gPar")
# dicrel<-dic[4]-dic
# 



