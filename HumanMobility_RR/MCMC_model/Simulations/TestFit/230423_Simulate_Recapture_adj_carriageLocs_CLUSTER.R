setwd("/data/pam/team284/sb62/scratch/Migration/SouthAfrica/mobility_model/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
iters=10000
scale=.06
boot=1
#########Load in data
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

extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999
#### check fit of cluster run simulation
for (chain in 1){
  load(paste0("./MCMC_model/Simulations/output/carriageLocs/sim.ans.",scale,".",iters,".",boot,"_adj2.RData"))
  load(paste0("./MCMC_model/Simulations/output/carriageLocs/dat.in2",scale,".",iters,".",boot,"_adj2.RData"))
  load(paste0("./MCMC_model/Simulations/output/carriageLocs/dat.in.all",scale,".",iters,".",boot,"_adj2.RData"))

plot(sim.ans[,1])
ans<-mcmc(sim.ans,start=iters/3,end=iters)
plot(ans[,1])
summary(ans)$statistics[,1]
1-rejectionRate(ans)
posteriors<-as.matrix(ans)
nsim<-4
pars<-matrix(nrow=nsim,ncol=ncol(ans))
for(i in 1:nsim){
  simno<-sample(nrow(posteriors),1)
  pars[i,]<-posteriors[simno,]
}
maxGen=30

true.modfit<-sub.modfit<-matrix(ncol=nsim,nrow=maxGen)
for(post in 1:nsim){
  # par<-pars[post,]

  par<-summary(ans)$statistics[,1]
  parHome<-par[1]
  extTranMatDat.tmp$pars$homeSus<-parHome
  Timehome=exp(parHome)/(1+exp(parHome))
  nInfecLoc<-rep(1,9)
  nInfecLoc[1:8]<-exp(par[2:9])
  nInfecLoc<-nInfecLoc/sum(nInfecLoc)
  ###No. detected
  avNoDetectByloc2<-colSums(rbind(hist(dat.in2$loc1,plot=FALSE,breaks=seq(0.5,9.5,1))$counts,hist(dat.in2$loc2,plot=FALSE,breaks=seq(0.5,9.5,1))$counts))
  # avNoDetectByloc2<-as.numeric(table(c(dat.in2$loc1,dat.in2$loc2)))
  no.detect<-avNoDetectByloc2
  probByloc<-no.detect/nInfecLoc
  probByloc<-probByloc/sum(probByloc)
  

  #####Truth
  nInfecLoc.true<-colSums(rbind(hist(dat.in.all$loc1,plot=FALSE,breaks=seq(0.5,9.5,1))$counts,hist(dat.in.all$loc2,plot=FALSE,breaks=seq(0.5,9.5,1))$counts))
  nInfecLoc.true<-nInfecLoc.true/sum(nInfecLoc.true)
  incidence.true<-nInfecLoc.true/extTranMatDat.tmp$popbyCell
  #####Check how it looks against population size
  incidence<-nInfecLoc/extTranMatDat.tmp$popbyCell
  tmp1<-cbind(exp(posteriors[,2:9]),1)
  tmp<-sweep(tmp1,1,rowSums(tmp1),"/")
  quantiles.probInfec<-apply(tmp,2,quantile,probs=c(0.025,0.5,0.975))
  
  
  ####Plot infection probability against the truth
  # pdf(file="/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/plots/infecProbLoc.pdf",height=2,width=2)
  # # quartz(width=2,height=2)
  # par(mar=c(3,3,1,1))
  plot(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[2,],pch=20,axes=F,ylim=c(0,0.5))
  points(extTranMatDat.tmp$popbyCell/1000000,nInfecLoc.true,col="red",pch=20,axes=F,ylim=c(0,0.5))
  arrows(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[1,],extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[3,],length=0)
  axis(1,cex=0.75,tck=-0.04,padj=-1.5,cex.axis=0.75)
  axis(2,cex=0.75,tck=-0.04,hadj=0.5,cex.axis=0.75,las=1)
  mtext(2,text="Proportion of infections",cex=0.75,line=1.7)
  mtext(1,text="Population size (x10^6)",cex=0.75,line=1.2)
  # dev.off()
  
  # }####end of ifelse for single par vs. multiple par
  ################################
  
  nlocs=nloc=9
  npop<-pop_2019
  pHome<-npop/sum(npop)
  min.range=-0.04
  max.range=0.6
  tmp.pHome.MUNIC<-pop2019.town
  # tmp.pHome.MUNIC<-tmp.pHome.MUNIC/sum(tmp.pHome.MUNIC)
  tmp.pHome<-extTranMatDat.tmp$popbyCell
  tmp.pHome<-tmp.pHome/sum(tmp.pHome)
  
  
  #########MUNICIPALITY
  ### Create mobility matrix for susceptible individuals
  tmpbase_pre<-cdr.mat.town
  tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
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
  move4<-sweep(move3,2,tmp.pHome,"*")
  TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")
  
  ##########
  probByGenA<-probByGenB<-rep(1,maxGen)
  gensB<-gensA<-(1:maxGen)
  maxgen.tmp<-max(c(gensA,gensB))
  TranMatArray<-array(NA,c(nlocs,nlocs,maxgen.tmp))
  TranMatArray[,,1]<-TranMat.tmp2
  if(Eigencpp){
    for (j in 2:maxgen.tmp){TranMatArray[,,j]<-eigenMapMatMult(TranMatArray[,,j-1], TranMat.tmp2)
    }
  }else{
    for (j in 2:maxgen.tmp){TranMatArray[,,j]<-TranMatArray[,,j-1]%*% TranMat.tmp2}
  }
  
  
  ### by province population size
  mrcaVec<-pHome
  
  ### true data counts for starting location.
  mrcaVecTRUE<-hist(dat.in.all$start,breaks=(0:nloc)+0.5,plot=F)$counts
  mrcaVecTRUE<-mrcaVecTRUE/sum(mrcaVecTRUE)
  # mrcaVec<-mrcaVecTRUE	
  
  ### destination infection probability
  TranMatArrayA.PreSamp<-TranMatArrayA<-TranMatArray[,,gensA]
  TranMatArrayA<-sweep(TranMatArrayA,2,probByloc,"*")####MULTIPLE PARAMETERS
  
  
  TranMatArrayB<-TranMatArray[,,gensB]
  TranMatArrayB<-sweep(TranMatArrayB,2,probByloc,"*")#####MULTIPLE PARAMETERS
  ###origin population
  TranMatArrayB.2<-sweep(TranMatArrayB,1,mrcaVec,"*")
  TranMatArray.PreSamp<-sweep(TranMatArrayA.PreSamp,1,mrcaVec,"*")
  TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
  
  ### sweep by true data sample:
  TranMatArrayA2.mrcaTRUE<-sweep(TranMatArrayA,1,mrcaVecTRUE,"*")

  
  TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
  TranMatArrayA3.mrcaTRUE<-matrix(TranMatArrayA2.mrcaTRUE,nlocs,nlocs*length(gensA))
  TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)
  probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)
  probAllPrs.mrcaTRUE<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3.mrcaTRUE)
  
  nogens<-rep(1:maxGen,each=nlocs)
  reflocs<-rep(1:nlocs,maxGen)
  nogens.mat<-outer(nogens,nogens,"+")
  a<-as.matrix(expand.grid(reflocs,reflocs))
  refDist<-pairwise_geodist[a]
  refDist.mat<-matrix(refDist,maxGen*nlocs,maxGen*nlocs)
  
  possGens<-expand.grid(1:maxGen,1:maxGen)
  possGenstot<-possGens[,1]+possGens[,2]
  
  distgenfromMRCA2<-distgenfromMRCA<-distgenfromMRCA.mrcaTRUE<-rep(NaN,maxGen)
  for (i in 1:maxGen){
    distgenfromMRCA[i]<-weighted.mean(pairwise_geodist,TranMatArrayA2[,,i])
    distgenfromMRCA.mrcaTRUE[i]<-weighted.mean(pairwise_geodist,TranMatArrayA2.mrcaTRUE[,,i])
    distgenfromMRCA2[i]<-weighted.mean(pairwise_geodist,TranMatArray.PreSamp[,,i])
  }
  
  weighted.distgen<-rep(NaN,maxGen)
  for (i in 1:maxGen){
    a<-which(nogens.mat==i)
    weighted.distgen[i]<-weighted.mean(refDist.mat[a],probAllPrs[a])
  }
  
  weighted.distgen2.mrcaTRUE<-weighted.distgen2<-rep(NaN,maxGen)
  for(i in 2:maxGen){
    a<-which(possGenstot==i)
    tmpWt2<-tmpWt<-rep(NaN,length(a))
    
    for (j in 1:length(a)){
      indA<-possGens[a[j],1]
      indB<-possGens[a[j],2]
      
      probDestA<-TranMatArray[,,indA]
      probDestB<-TranMatArray[,,indB]
      
      
      probDestA<-sweep(probDestA,2,probByloc,"*") ##MULTIPLE PARS
      probDestB<-sweep(probDestB,2,probByloc,"*")
      
      
      dist.tmp<-rep(NaN,9)
      wtDist2<-wtDist<-matrix(NaN,0,2)
      for(k in 1:9){
        b<-outer(probDestA[k,],probDestB[k,],"*")*mrcaVec[k]
        # b<-b/sum(b)
        # dist.tmp[k]<-weighted.mean(pairwise_geodist,b)
        wtDist<-rbind(wtDist,cbind(as.vector(pairwise_geodist),as.vector(b)))
        b<-outer(probDestA[k,],probDestB[k,],"*")*mrcaVecTRUE[k]
        wtDist2<-rbind(wtDist2,cbind(as.vector(pairwise_geodist),as.vector(b)))
      }
      # tmpWt[j]<-mean(dist.tmp)
      tmpWt[j]<-weighted.mean(wtDist[,1],wtDist[,2])
      tmpWt2[j]<-weighted.mean(wtDist2[,1],wtDist2[,2])
      # tmpWt[j]<-weighted.mean(dist.tmp,mrcaVec)
    }
    # weighted.distgen2[i]<-weighted.mean(wtDist[,1],wtDist[,2])
    weighted.distgen2[i]<-mean(tmpWt)
    weighted.distgen2.mrcaTRUE[i]<-mean(tmpWt2)
  }
  sub.modfit[,post]<-weighted.distgen2
  true.modfit[,post]<-distgenfromMRCA2
  print(post)
}
## Distance#####################
dat.in.tmp<-dat.in2
GenSeqMax<-seq(1,50,1)
GenSeqMin<-GenSeqMax-2
GenSeqMin[which(GenSeqMin<0)]<-0
GenSeqMid<-(GenSeqMin+GenSeqMax)/2

distByGenStart2<-distByGenStart1<-distByGenSamp<-distByGen<-rep(NaN,length(GenSeqMax))
for(i in 1:length(GenSeqMin)){
  a<-which(dat.in.all$meanGen>=GenSeqMin[i]&dat.in.all$meanGen<GenSeqMax[i])
  distByGen[i]<-mean(dat.in.all$dist[a])
  
  a<-which(dat.in.all$nogens>=GenSeqMin[i]&dat.in.all$nogens<GenSeqMax[i])
  distByGenStart1[i]<-mean(dat.in.all$distStart1[a])
  
  a<-which(dat.in.all$nogens2>=GenSeqMin[i]&dat.in.all$nogens2<GenSeqMax[i])
  distByGenStart2[i]<-mean(dat.in.all$distStart2[a])
  
  a<-which(dat.in.tmp$meanGen>=GenSeqMin[i]&dat.in.tmp$meanGen<GenSeqMax[i])
  distByGenSamp[i]<-mean(dat.in.tmp$dist[a])
}

data.sim <- cbind(GenSeqMid,distByGen,distByGenSamp,distByGenStart1,distByGenStart2)
data.sim <- data.table(data.sim)
colnames(data.sim) <- c("GenSeqMid","distByGen","distByGenSamp","distByGenStart1","distByGenStart2")

modfit <- cbind(c(1:maxGen),t(apply(sub.modfit, 1, function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm = TRUE))), 
                t(apply(true.modfit, 1, function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm = TRUE))))
# modfit <- cbind(c(1:maxGen),distgenfromMRCA2,weighted.distgen,weighted.distgen2,weighted.distgen2.mrcaTRUE)
# colnames(modfit) <- c("gens","distgenfromMRCA2","weighted.distgen","weighted.distgen2","weighted.distgen2.mrcaTRUE")

modfit <- data.table(modfit)
colnames(modfit) <- c("gens","weighted.distgen2.lower","weighted.distgen2","weighted.distgen2.upper",
                      "distgenfromMRCA2.lower","distgenfromMRCA2","distgenfromMRCA2.upper")

simFit <-ggplot()+
  geom_point(data=data.sim,aes(x=GenSeqMid,y=distByGen,color="black"),color="black",alpha=0.7) +
  geom_point(data=data.sim,aes(x=GenSeqMid,y=distByGenSamp,color="red"),color="red",alpha=0.7)+
  
  geom_line(data=modfit,aes(x=gens,y=distgenfromMRCA2,color="black"),lwd=1.5,color="black",linetype="dashed")+
  geom_ribbon(data=modfit,aes(ymin=weighted.distgen2.lower,ymax=weighted.distgen2.upper,x=gens),alpha=0.2,fill="red")+
  geom_line(data=modfit,aes(x=gens,y=weighted.distgen2,color="red"),lwd=1.5,color="red",linetype="dashed")+
  geom_ribbon(data=modfit,aes(ymin=distgenfromMRCA2.lower,ymax=distgenfromMRCA2.upper,x=gens),alpha=0.4,fill="grey")+
  
  # geom_line(data=modfit,aes(x=gens,y=weighted.distgen2.mrcaTRUE,color="blue"),lwd=1.5,color="blue",linetype="dashed")+
  # geom_line(data=modfit,aes(x=gens,y=weighted.distgen,color="orange"),lwd=1.5,color="orange",linetype="dashed")+
  
  labs(color="Type")+
  theme_bw()+
  scale_x_continuous(limits = c(0, 50), name = "Evolutionary Time (generations)",    # Features of the first axis
                     sec.axis = sec_axis( trans=~.*35/365, name="Evolutionary Time (years)") )+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  ylab("Spatial Distance (km)")
simFit
# ggsave(filename = paste0("./MCMC_model/Simulations/Plots/simFit",chain,".pdf"))
}


simFit
