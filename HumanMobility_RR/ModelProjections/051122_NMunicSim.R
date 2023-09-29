#### Simulating epidemics from each municipality

library(data.table)
library(ggplot2)
library(dplyr)
load("./modelinput_data/cdr.mat.town.one.RData")
cdr.mat.town<-cdr.mat.town.one
load("./modelinput_data/pairwise_geodist.town.RData") 
load("./modelinput_data/pop_municipality.2017LS.RData") 
## Functions
SeasonFunc<-function(day,r0=1.1,amp=0.03){
  d<-r0+amp*(sin(-pi/2+2*pi*day/365))
  return(d)}
mean=35
var=35^2
shape=mean^2/var
scale=var/mean
tranmat.cdr <- cdr.mat.town
##INPUT TRUE MOBILITY 
load("./ModelProjections/data/TranMatArray.234x234_adj.RData")
TranMatArray.1<-TranMatArray.234x234

provs <- c("Eastern Cape","Free State","Gauteng","KwaZulu-Natal","Limpopo","Mpumalanga","North West","Northern Cape","Western Cape")
tn <- rownames(cdr.mat.town)
pairwise_geodist.town <- pairwise_geodist.town[tn,tn]
pop2019.town <- pop2019.town[tn]
### Overall####
##Overall Function 
simFunction<-function(tranmat,specificStart=specificStart,noGens,overdis=F,R0=1,amp=0.15){
  nloc<-234
  if(length(specificStart)==0){
    startloc<-sample(nloc,1)}else{startloc=specificStart}
  genNo<-1
  whereFinal<-where<-startloc
  timeFinal<-time<-1
  nseed=1
  
  for (jj in 2:noGens){
    Reff<-SeasonFunc(time,r0=R0,amp=amp) ## With Seasonality normal
    # Reff=R0 ### Without Seasonality
    nOffspring<-rpois(nseed,Reff)
    # nOffspring<-Reff
    if(overdis){nOffspring<-rnbinom(nseed,mu=Reff,size=1)}
    tmp2<-tmp<-NULL
    if(sum(nOffspring)==0)break
    for (ll in 1:nseed){
      tmp<-c(tmp,sample(nloc,nOffspring[ll],replace=T,prob=tranmat[where[ll],]))
      tmp2<-c(tmp2,time[ll]+round(rgamma(nOffspring[ll],shape=shape,scale=scale)))
    }
    if (length(tmp2) > 30000) {print("SKIP (too long)"); next}
    # if (length(tmp2) > 50000) {print("SKIP (too long)"); next}
    
    where<-tmp
    time<-tmp2
    nseed=length(where)
    whereFinal<-c(whereFinal,where)
    timeFinal<-c(timeFinal,time)
    genNo<-c(genNo,rep(jj,length(where)))
  }
  
  out.dat<-data.frame(loc=whereFinal,time=timeFinal,gen=genNo)
  dim(out.dat)
  return(out.dat)
}
noGens=120
## time in days, generation of that day.
nboot=200
mat.boots<-list()
plot.nmuic.list <- list()
plot.dist.list <- list()
# for (boot in 1:234){
for (boot in 1:nboot) {
  times <- NULL
  repeat {
    # specificStart <- boot
    specificStart <- sample(234,1,prob=c(pop2019.town))
    sim<-simFunction(tranmat=tranmat.cdr,specificStart=specificStart,noGens=noGens,overdis=F,R0=1,amp=0.15) 
    print(dim(sim)[1]) 
    a <- sim$time[length(sim$time)]
    times <- c(times,a)
    if ( sim$time[length(sim$time)]>4380  ){break} ## Overall
    # if ( sim$time[length(sim$time)]>2000  ){break}  ## NVT 2010 after PCV
  }
  
  
  if(dim(sim)[1]>50000){next}
  
  hist(sim$time,breaks=c(100))
  simOut<-sim
  if ( max(simOut$time) > 18250) { ngen = 18250 } else( ngen=max(simOut$time))
  npairs<-which(simOut$time==ngen)[1]
  simOut.tmp2 <- simOut[1:npairs,]
  days <- sort(unique(simOut.tmp2$time))
  ts <- hist(days,c(30),plot=F)$breaks
  dmin <- ts[1:(length(ts)-1)]
  dmax <-ts[2:length(ts)]
  print(length(dmax))
  mat.tmp <- matrix(nrow=npairs,ncol=8)
  ## Weigthed mean distance across all pairs 
  for(i in 1:npairs ){
    mat.tmp[i,1] <- weighted.mean(pairwise_geodist.town[simOut.tmp2[i,1],], TranMatArray.1[simOut.tmp2[i,1],,simOut.tmp2[i,3]] )
    mat.tmp[i,2] <- simOut.tmp2[i,3]
    mat.tmp[i,3] <- simOut.tmp2[i,2]
    mat.tmp[i,4] <- tn[simOut.tmp2[i,1]]
    mat.tmp[i,5] <-pop2019.town[simOut.tmp2[i,1]]
    ## N Munic per Gen
    tmp.gen <- simOut.tmp2[1:i,]
    mat.tmp[i,6] <- length(table(unique(tmp.gen$loc)))
    ## proportion home 
    start <- simOut.tmp2$loc[1]
    mat.tmp[i,7] <- length(which(simOut.tmp2[1:i,]$loc==start))/i
    mat.tmp[i,8] <- length(which(simOut.tmp2[1:i,]$loc==20))/i 
  }
  mat.tmp <- data.table(mat.tmp)
  colnames(mat.tmp) <-c("Distance","gens","days","municipality","population","NperGen","propHome","propJoBurg")
  mat.tmp$Distance<-as.numeric(mat.tmp$Distance)
  mat.tmp$gens<-as.numeric(mat.tmp$gens)
  mat.tmp$days<-as.numeric(mat.tmp$days)
  mat.tmp$population<-as.numeric(mat.tmp$population)
  mat.tmp$NperGen <- as.numeric(mat.tmp$NperGen)
  mat.tmp$propHome <- as.numeric(mat.tmp$propHome)
  mat.tmp$propJoBurg <- as.numeric(mat.tmp$propJoBurg)
  
  
  ### Mean at each generation
  maxGen=max(mat.tmp$gens,na.rm = T)
  mat <-matrix(nrow=maxGen,ncol=7)
  for ( j in 1:maxGen){
    tmp <- subset(mat.tmp, mat.tmp$gens==j )
    mat[j,1] <- mean(tmp$Distance)
    # mat[j,1] <- median(tmp$Distance)
    
    mat[j,2] <- mean(tmp$days)
    mat[j,3] <- mat.tmp$municipality[1]
    mat[j,4] <- mat.tmp$population[1]
    mat[j,5] <-  mean(tmp$NperGen)
    # mat[j,5] <-  median(tmp$NperGen)
    
    mat[j,6] <-  mean(tmp$propHome)
    mat[j,7] <-  mean(tmp$propJoBurg)
    
  }
  
  mat <- data.table(mat)
  colnames(mat) <-c("Distance","days","municipality","population","NperGen","propHome","propJoBurg")
  mat$gens <- c(1:maxGen)
  mat$Distance<-as.numeric(mat$Distance)
  mat$gens<-as.numeric(mat$gens)
  mat$days<-as.numeric(mat$days)
  mat$population<-as.numeric(mat$population)
  mat$NperGen<-as.numeric(mat$NperGen)
  mat$propHome<-as.numeric(mat$propHome)
  mat$propJoBurg<-as.numeric(mat$propJoBurg)
  
  mat$boot <- boot
  mat.boots[[boot]] <- mat
  print(paste0("We are at",boot,"iteration"))
}
# mat.boots<-mat.boots[1:234]
mat.tot <- rbindlist(mat.boots)

### mean at each generation
mat.tot.fin <- matrix(nrow=noGens,ncol=16)
for (i in 1:noGens){
  tmp <- subset(mat.tot, mat.tot$gens==i)
  mat.tot.fin[i,1] <- mean(tmp$NperGen)
  mn <- mean(tmp$NperGen)
  a <- sd(tmp$NperGen)
  b <- sd(tmp$NperGen)/sqrt(nrow(tmp))
  mat.tot.fin[i,2] <- mn-(1.96*b)
  mat.tot.fin[i,3] <- mn+(1.96*b)
  mat.tot.fin[i,1] <- mean(tmp$NperGen)
  
  
  mat.tot.fin[i,4] <- mean(tmp$Distance)
  mn <- mean(tmp$Distance)
  a <- sd(tmp$Distance)
  b <- sd(tmp$Distance)/sqrt(nrow(tmp))
  mat.tot.fin[i,5] <- mn-(1.96*b)
  mat.tot.fin[i,6] <- mn+(1.96*b)
  
  mat.tot.fin[i,7] <- mean(tmp$propHome)
  mn <- mean(tmp$propHome)
  a <- sd(tmp$propHome)
  b <- sd(tmp$propHome)/sqrt(nrow(tmp))
  mat.tot.fin[i,8] <- mn-(1.96*b)
  mat.tot.fin[i,9] <- mn+(1.96*b)
  
  mat.tot.fin[i,10] <- mean(tmp$propJoBurg)
  mn <- mean(tmp$propJoBurg)
  a <- sd(tmp$propJoBurg)
  b <- sd(tmp$propJoBurg)/sqrt(nrow(tmp))
  mat.tot.fin[i,11] <- mn-(1.96*b)
  mat.tot.fin[i,12] <- mn+(1.96*b)
  mat.tot.fin[i,13] <- i
  mat.tot.fin[i,14] <- mean(tmp$days)
  mat.tot.fin[i,15] <- median(tmp$NperGen)
  mat.tot.fin[i,16] <- median(tmp$Distance)
  
  
  
}
mat.tot.fin.overall <-  data.table(mat.tot.fin)
colnames(mat.tot.fin.overall) <- c("NperGen","NperGen_lower","NperGen_upper",
                                   "Distance","Distance_lower","Distance_upper",
                                   "propHome","propHome_lower","propHome_upper",
                                   "propJoburg","propJoburg_lower","propJoburg_upper","gen","days","medNperGen","medDistance")

save(mat.tot.fin.overall,file="./ModelProjections/data/mat.tot.fin.overall_adj.RData")
save(mat.tot,file="./ModelProjections/data/mat.tot_adj.RData")


######PLOTS
plot.dist <- ggplot() +
  geom_line( data=mat.tot,aes(days/365, Distance,group=as.character(boot)),color="grey")+
  geom_line(data=mat.tot.fin.overall,aes(days/365, medDistance),color="black" )+
  theme_bw()+
  xlab("Years")+
  ylab("Distance (km)")+
  ggtitle("A")+
  theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))
plot.nmuic <- ggplot() +
  geom_line(data=mat.tot,aes(days/365, NperGen,group=as.character(boot)) ,color="grey")+
  geom_line(data=mat.tot.fin.overall,aes(days/365, medNperGen),color="black" )+
  theme_bw()+
  scale_color_discrete(name="iter")+
  theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
  ylim(1,10)+
  ggtitle("B")+
  
  xlab("Years")+
  ylab("Municipalities\n Visited (N)")
plot.probs.overall <- ggplot() +
  geom_line(data=mat.tot,aes(days/365, propHome,group=as.character(boot)),color="grey")+
  geom_line(data=mat.tot.fin.overall,aes(days/365, propHome),color="black")+
  theme_bw()+
  theme(axis.text = element_text(size=20),axis.title = element_text(size=20))+
  xlab("Years")+
  ggtitle("C")+
  
  ylab("Proportion in \nStart Municipality")

library(patchwork)

(plot.dist + plot.nmuic )#+ plot.probs.overall)

# ggsave("./ModelProjections/plots/sims.overall.dist.Nmunic.pdf",width=12,height=4)


