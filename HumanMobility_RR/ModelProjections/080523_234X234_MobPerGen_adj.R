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
load('./modelinput_data/dat.tmp.allser.RData')
load('./modelinput_data/pairwise_geodist.RData')
load("./modelinput_data/pop_municipality.2017LS.RData") 
load("./modelinput_data/pairwise_geodist.town.RData") 
load('./modelinput_data/cdr.mat.town.one.RData')
cdr.mat.town<-cdr.mat.town.one

maxGen=250
## set provs order
provs<-colnames(cdr.mat)
nloc.munic<-length(pop2019.town)
nam.b<-colnames(cdr.mat.town)
pop2019.town<-pop2019.town[nam.b]
pairwise_geodist.town<-pairwise_geodist.town[nam.b,nam.b]

#### Index municipalities and provinces
splitnames <- matrix(nrow=nloc.munic,ncol=2)
for (us in 1:nloc.munic) {
  for (prov in 1:2) { 
    splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
  }
}

### Collapse province populations by municipality populations
pop_2019<-vector(mode="numeric",length=9)
for (sp in 1:length(provs)){
  sel_sp <- which(splitnames[,2]==provs[sp])
  pop_2019[sp]<-sum(pop2019.town[sel_sp])
}
names(pop_2019) <- provs

nam.a<-colnames(cdr.mat)
pop_2019<-pop_2019[nam.a]
pairwise_geodist<-pairwise_geodist[nam.a,nam.a]


nlocs=9
Eigencpp=TRUE
calcAllProbs=FALSE

# maxGen=maxTranGens
sourceCpp("./MCMC_model/MatrixMultiplication.cpp")
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic_adj.R")

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
### all chains

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

## true detect
NoDetectByYearGPSC<-table(c(dat.inMaster$loc1,dat.inMaster$loc2),c(dat.inMaster$GPSC,dat.inMaster$GPSC),c(dat.inMaster$year1,dat.inMaster$year2))
avNoDetectByloc2<-as.numeric(table(c(dat.in2$loc1,dat.in2$loc2)))
no.detect<-avNoDetectByloc2

##### Adjust to account for year variability
no.detect<-apply(rbind(apply(table(dat.in2$YearREF1,dat.in2$loc1), 2, mean), apply(table(dat.in2$YearREF2,dat.in2$loc2), 2, mean)), 2, sum)
fit.par<-apply(rbind(summary(chain1)$statistics[,1],summary(chain2)$statistics[,1],summary(chain3)$statistics[,1]),2,mean)
# fit.par<-summary(ans)$statistics[,1]
fit.par[1]
nInfecLoc<-rep(1,9)
nInfecLoc[1:8]<-exp(fit.par[2:9])
nInfecLoc<-nInfecLoc/sum(nInfecLoc)
probByloc<-no.detect/nInfecLoc


nloc=nlocs=9
tmp.pHome<-pop2019.town
# tmp.pHome.PROV<-extTranMatDat.tmp$popbyCell

# avNoDetectByloc2<-as.numeric(table(c(dat.in2$loc1,dat.in2$loc2)))

extTranMatDat.tmp$pars$homeSus<-fit.par[1]

# ### Create mobility matrix for susceptible individuals
# tmpbase<-cdr.mat.town
# tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
# tmppar <- min.range+tmppar1*(max.range-min.range)
# tmpdiag<-diag(tmpbase)-tmppar
# tmpdiag[which(tmpdiag>0.99999)]<-0.99999
# diag(tmpbase)<-0
# tmpbase_pre2<-sweep(tmpbase_pre2,1,rowSums(tmpbase_pre2),"/")
# tmpbase_pre2<-sweep(tmpbase_pre2,1,(1-tmpdiag)/(1-diag(tmpbase_pre2)),"*")
# diag(tmpbase_pre2)<-tmpdiag
# weighted.mean(tmpdiag,pop2019.town)
# mean(tmpdiag)
# nloc.munic=nrow(tmpbase_pre2)

### Create mobility matrix for susceptible individuals
tmpbase_pre<-cdr.mat.town
tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
tmppar <- min.range+tmppar1*(max.range-min.range)
tmpdiag<-diag(tmpbase_pre)-tmppar
tmpdiag[which(tmpdiag>0.99999)]<-0.99999
diag(tmpbase_pre)<-0
tmpbase_pre1<-sweep(tmpbase_pre,1,rowSums(tmpbase_pre),"/")
tmpbase_pre2<-sweep(tmpbase_pre1,1,(1-tmpdiag)/(1-diag(tmpbase_pre1)),"*")
mean(tmpdiag)
diag(tmpbase_pre2)<-tmpdiag
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

# setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/")
# save(tmpbase,file="munic.adj.mob.RData")
# 
move3<-tcrossprod(tmpbase,tmpbase)
move4<-sweep(move3,2,tmp.pHome,"*")
TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")

# move3<-tcrossprod(tmpbase,tmpbase)
# move4<-sweep(move3,2,tmp.pHome,"*")
# move5<-sweep(move4,1,tmp.pHome,"*")
# TranMat.tmp2<-sweep(move5,1,rowSums(move5),"/")


mean(diag(TranMat.tmp2))
# weighted.mean(diag(TranMat.tmp2),pop2019.town)

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
TranMatArray.234x234<-TranMatArray.1
# save(TranMatArray.234x234,file="./ModelProjections/data/TranMatArray.234x234_adj.RData")

tmpbase.1 <- apply(tmpbase, 1, function(x){
  split(x,splitnames[,2]) %>% sapply(sum)
})  %>% transpose()

tmpbase2 <- apply(tmpbase.1,2,function(x) {
  mapply( weighted.mean
          , x = split(x,splitnames[,2])
          , w = split(pop2019.town,splitnames[,2])
  )
})
mean(diag(tmpbase2))
# weighted.mean(diag(tmpbase2),pop_2019)

mean(diag(tmpbase))
# weighted.mean(diag(tmpbase),pop2019.town)


tmpbase.1 <- apply(TranMat.tmp2, 1, function(x){
  split(x,splitnames[,2]) %>% sapply(sum)
})  %>% transpose()

tmpbase2 <- apply(tmpbase.1,2,function(x) {
  mapply( weighted.mean
          , x = split(x,splitnames[,2])
          , w = split(pop2019.town,splitnames[,2])
  )
})

onegendiag<-diag(TranMatArray.234x234[,,1])
names(onegendiag)<-names(pop2019.town)
