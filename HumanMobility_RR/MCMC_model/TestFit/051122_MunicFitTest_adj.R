#######Test Fit of parameters when fit at municipality level##########
###### Adjusting down to 9X9 BEFORE MATRIX MULTIPLICATION --- MATCHING THE WAY WE FIT
##### Adding variables from other scripts and loading libraries
iters=20000
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
### Adding variables
# maxGen=600
# maxTranGens=maxGen

nlocs=9
Eigencpp=TRUE
calcAllProbs=FALSE

# maxGen=maxTranGens
sourceCpp("./MCMC_model/MatrixMultiplication.cpp")
source("./MCMC_model/LikelihoodFunctions/LikFunc.mcmc.Munic.R")


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

###scale 0.08

### all chains
load(paste0("./MCMC_model/outputs/mechan_model/ans.munic1.",iters,".08_adj.RData"))
chain1<-ans.munic
chain1<-chain1[1000:4000,]
load(paste0("./MCMC_model/outputs/mechan_model/ans.munic2.",iters,".08_adj.RData"))
chain2<-ans.munic
chain2<-chain2[1000:4000,]
load(paste0("./MCMC_model/outputs/mechan_model/ans.munic3.",iters,".08_adj.RData"))
chain3<-ans.munic
chain3<-chain3[1000:4000,]
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
# no.detect<-apply(rbind(apply(table(dat.in2$YearREF1,dat.in2$loc1), 2, median), apply(table(dat.in2$YearREF2,dat.in2$loc2), 2, median)), 2, sum)
no.detect<-apply(rbind(apply(table(dat.in2$YearREF1,dat.in2$loc1), 2, mean), apply(table(dat.in2$YearREF2,dat.in2$loc2), 2, mean)), 2, sum)

### input pars into fit.par to plot uncertainty around model estimate
maxGen=500

# maxGen=1000
nboot=10
true.dat<-dat<-matrix(nrow=maxGen,ncol=ncol(pars))
for (boot in 1:ncol(pars)){
  fit.par<-pars[boot,]
  fit.par<-apply(ans.munic.chains,2,mean)
fit.par[1]
nInfecLoc<-rep(1,9)
nInfecLoc[1:8]<-exp(fit.par[2:9])
nInfecLoc<-nInfecLoc/sum(nInfecLoc)
probByloc<-no.detect/nInfecLoc
# probByloc<-probByloc/sum(probByloc)
incidence<-nInfecLoc/pop_2019

# 
tmp1<-cbind(exp(posteriors[,2:9]),1)
tmp<-sweep(tmp1,1,rowSums(tmp1),"/")
quantiles.probInfec<-apply(tmp,2,quantile,probs=c(0.025,0.5,0.975))
df<-rbind(quantiles.probInfec,extTranMatDat.tmp$popbyCell/1000000)

# # Plot the incidence vs. the populations size to test detection parameters
# par(mar=c(3,3,1,1))
# plot(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[2,],pch=20,axes=F,ylim=c(0,0.5))
# arrows(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[1,],extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[3,],length=0)
# axis(1,cex=0.75,tck=-0.04,padj=-1.5,cex.axis=0.75)
# axis(2,cex=0.75,tck=-0.04,hadj=0.5,cex.axis=0.75,las=1)
# mtext(2,text="Proportion of infections",cex=0.75,line=1.7)
# mtext(1,text="Population size (x10^6)",cex=0.75,line=1.2)
df.propperprov<-df

df.propperprov<-data.table(t(df.propperprov))
colnames(df.propperprov)<-c("lowerprop","prop","upperprop","popsize")
ggplot(df.propperprov,aes(x=popsize,y=prop))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lowerprop,ymax=upperprop))+
  theme_classic()+
  xlab("Population size (x10^6)")+
  ylab("Proportion of Infections")+
  theme(axis.text = element_text(size=20),axis.title=element_text(size=20))


# summary(lm(df$prop~df$popsize))
# ggsave("./MCMC_model/TestFit/Plots/propperpop_adj.pdf",width=4,height=4)


nloc=nlocs=9
# tmp.pHome<-pop2019.town
tmp.pHome.PROV<-extTranMatDat.tmp$popbyCell

# avNoDetectByloc2<-as.numeric(table(c(dat.in2$loc1,dat.in2$loc2)))

### Set maxGen for plot
extTranMatDat.tmp$pars$homeSus<-fit.par[1]

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
weighted.mean(diag(tmpbase2),pop_2019)

move3<-tcrossprod(tmpbase2,tmpbase2)
move4<-sweep(move3,2,tmp.pHome.PROV,"*")
TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")

mean(diag(TranMat.tmp2))
weighted.mean(diag(TranMat.tmp2),pop_2019)

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

# save(TranMatArray,file="./MCMC_model/TestFit/Data/TranMatArray.munic.RData")

quantile(diag(TranMatArray[,,12]),probs=c(0.025,0.5,0.975))
# no.detect<-avNoDetectByloc2

# nInfecLoc<-rep(1,9)
# nInfecLoc[1:8]<-exp(fit.par[2:9])
# probByloc<-no.detect/nInfecLoc
# probByloc<-probByloc/sum(probByloc)
# mrcaVec2<-mrcaVec<-extTranMatDat.tmp$popbyCell
mrcaVec2<-mrcaVec<-pop_2019


TranMatArrayA<-TranMatArray[,,gensA]
TranMatArrayA.BASE<-TranMatArrayA<-sweep(TranMatArrayA,2,probByloc,"*")

TranMatArrayB.BASE<-TranMatArrayB<-TranMatArray[,,gensB]
TranMatArrayB<-sweep(TranMatArrayB,2,probByloc,"*")

TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
# TranMatArrayB.2<-sweep(TranMatArrayB,1,mrcaVec,"*")
TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)

probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)


nogens<-rep(1:maxGen,each=nlocs)
reflocs<-rep(1:nlocs,maxGen)
nogens.mat<-outer(nogens,nogens,"+")
a<-as.matrix(expand.grid(reflocs,reflocs))
refDist<-pairwise_geodist[a]
refDist.mat<-matrix(refDist,maxGen*nlocs,maxGen*nlocs)

possGens<-expand.grid(1:maxGen,1:maxGen)
possGenstot<-possGens[,1]+possGens[,2]

distgenfromMRCA2<-distgenfromMRCA<-rep(NaN,maxGen)
for (i in 1:maxGen){
  distgenfromMRCA[i]<-weighted.mean(pairwise_geodist,TranMatArrayA.BASE[,,i])
  distgenfromMRCA2[i]<-weighted.mean(pairwise_geodist,TranMatArrayB.BASE[,,i])
}

weighted.distgen<-rep(NaN,maxGen)
for (i in 1:maxGen){
  a<-which(nogens.mat==i)
  weighted.distgen[i]<-weighted.mean(refDist.mat[a],probAllPrs[a])
}

weighted.distgen[1]<-0
dat[,boot]<-weighted.distgen
true.dat[,boot]<-distgenfromMRCA
print(boot)
}
dat.quants<-matrix(nrow=maxGen,ncol=4)
dat.quants[,1]<-apply(dat,1,mean)
dat.quants[,2:4]<-t(apply(dat,1,function(x) quantile(x,probs=c(0.025,0.5,0.975))))
dat.quants<-data.table(dat.quants)
colnames(dat.quants)<-c("mean","lowerCI","median","upperCI")
dat.quants$gens<-1:maxGen
dat.quants$years<-(dat.quants$gens*35)/365
save(dat.quants,file="./MCMC_model/TestFit/Data/dat.quants.Munic.RData")
##### ignoring sampling probability
true.dat.quants<-matrix(nrow=maxGen,ncol=4)
true.dat.quants[,1]<-apply(true.dat,1,mean)
true.dat.quants[,2:4]<-t(apply(true.dat,1,function(x) quantile(x,probs=c(0.025,0.5,0.975))))
true.dat.quants<-data.table(true.dat.quants)
colnames(true.dat.quants)<-c("mean","lowerCI","median","upperCI")
true.dat.quants$gens<-1:maxGen
true.dat.quants$years<-(true.dat.quants$gens*35)/365
# save(true.dat.quants,file="./MCMC_model/TestFit/Data/true.dat.quants.Munic.RData")


# ################################################################################################
load("./MCMC_model/TestFit/Data/overall_data.RData")
# load("./MCMC_model/TestFit/Data/overall_data_multiTree.RData")


load("./MCMC_model/TestFit/Data/true.dat.quants.Munic_adj.RData")
load("./MCMC_model/TestFit/Data/dat.quants.Munic_adj.RData")

################################################################################################
# quartz(width=2,height=3.5)
# par(mar=c(3,3,1,1))
# plot(gDistMid,meanDist.emp.bs.qt[2,],ylim=c(0,500),xlim=c(0,100),type="n",axes=F)
# polygon(xpoly,ypoly1,col=rgb(t(col2rgb("light blue")/255),alpha=0.5),border=NA)
# lines(gDistMid,meanDist.emp.bs.qt[2,],col="blue")
# # lines(gDistMid,cummeanDist.emp.bs.qt[2,],col="blue",lty=2)
# lines((1:maxGen)*(2/12),weighted.distgen,col="red",lty=2)
# # lines((1:maxGen)*(2/12),dat.quants[,c("median")],col="red",lty=2)
# lines((1:maxGen)*(2/12),distgenfromMRCA,col="darkgreen",lty=2)
# axis(1,cex=0.75,tck=-0.04,padj=-1.5,cex.axis=0.75)
# axis(2,cex=0.75,tck=-0.04,hadj=0.5,cex.axis=0.75,las=1)
# mtext(2,text="Spatial distance (km)",cex=0.75,line=1.7)
# mtext(1,text="Evolution time (years)",cex=0.75,line=1.2)

##########Visualize Fit####
colors=c("Data"="blue","Model"="red","Implied Truth"="#BF40BF")
plot.munic<-ggplot()+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=dat.quants,aes(x=years,y=median,color="Model"))+
  geom_line(data=true.dat.quants,aes(x=years,y=median,color="Implied Truth"),linetype="dashed")+
  geom_line(data=overall_data,aes(x=years,y=median,color="Data"))+
  theme(axis.title = element_text(size=25),axis.text = element_text(size=20),
        legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
  geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.4, fill="blue")+
  geom_ribbon(data=dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.2, fill="red")+
  geom_ribbon(data=true.dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="#BF40BF")+
  scale_color_manual(values = colors,breaks=c("Data","Model","Implied Truth"),limits=c("Data","Model","Implied Truth"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  xlab("Evolutionary Time (years)")+
  theme(legend.position = "none")+
  
  ylab("Distance (km)")#+
  # xlim(0,3)+
  # ylim(0,300)
plot.munic
 
####JUST BLUE LINES
  ggplot()+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    # geom_line(data=dat.quants,aes(x=years,y=median,color="Model"))+
    # geom_line(data=true.dat.quants,aes(x=years,y=median,color="Implied Truth"))+
    geom_line(data=overall_data,aes(x=years,y=median,color="Data"))+
    theme(axis.title = element_text(size=25),axis.text = element_text(size=20),
          legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
    geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="blue")+
    # geom_ribbon(data=dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="red")+
    # geom_ribbon(data=true.dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="green")+
    scale_color_manual(values = colors,breaks=c("Data","Model","Implied Truth"),limits=c("Data","Model","Implied Truth"))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme(legend.position = "none")+
    xlab("Evolution Time (years)")+
    ylab("Distance (km)")#+
ggsave("./MCMC_model/TestFit/munic.RAWDATA.pdf",width=5,height=5.5)


####JUST BLUE LINES + Red
ggplot()+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=dat.quants,aes(x=years,y=median,color="Model"))+
  # geom_line(data=true.dat.quants,aes(x=years,y=median,color="Implied Truth"))+
  geom_line(data=overall_data,aes(x=years,y=median,color="Data"))+
  theme(axis.title = element_text(size=25),axis.text = element_text(size=20),
        legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
  geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="blue")+
  geom_ribbon(data=dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="red")+
  # geom_ribbon(data=true.dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="green")+
  scale_color_manual(values = colors,breaks=c("Data","Model","Implied Truth"),limits=c("Data","Model","Implied Truth"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme(legend.position = "none")+
  xlab("Evolution Time (years)")+
  ylab("Distance (km)")#+
ggsave("./MCMC_model/TestFit/Plots/munic.JUSTMODEL.pdf",width=5,height=5.5)


####Inset
plot.munic.inset<-ggplot()+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=dat.quants,aes(x=years,y=median,color="Model"))+
  geom_line(data=true.dat.quants,aes(x=years,y=median,color="Implied Truth"),linetype="dashed")+
  geom_point(data=overall_data,aes(x=years,y=median,color="Data"))+
  theme(axis.title = element_text(size=25),axis.text = element_text(size=20),
        legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
  # geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.4, fill="blue")+
  geom_errorbar(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.7, color="blue")+
  
  geom_ribbon(data=dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.2, fill="red")+
  # geom_ribbon(data=true.dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="green")+
  scale_color_manual(values = colors,breaks=c("Data","Model","Implied Truth"),limits=c("Data","Model","Implied Truth"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  xlab("Evolutionary Time (years)")+
  theme(legend.position = "none")+
  
  ylab("Distance (km)")+
xlim(0,3)+
  ylim(0,300)
plot.munic.inset
ggsave("./MCMC_model/TestFit/Plots/municfit.inset.pdf",width=5,height=5.5)
