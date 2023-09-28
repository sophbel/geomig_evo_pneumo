library(data.table)
library(ggplot2)
##########ORGANIZE THE DATA #####################

maxCarriageDuration<-365
sites<-c("Gambia","Kilifi")
clear.mat<-matrix(nrow=length(sites),ncol=7)
colnames(clear.mat)<-c("site","mean_dur","gentime_mean","sd","gentime_med","lowerCI","upperCI")
clear.mat[,1]<-sites

#########GAMBIA, CHRISPIN CHAGUZA FRONTIERS ########
nboot=10000
carr.mat.gamb<-matrix(nrow=nboot,ncol=2)
for(k in 1:nboot){
  durs<-rexp(1000,rate=0.0263)
  durs[which(durs>maxCarriageDuration)]<-NA
  durs<-durs[!is.na(durs)]
  
  d<-sample(durs,1,replace = FALSE,prob=durs)
  p<-runif(1,1,5)####incubation period
  t<-runif(1,0,p+d)
  # t<-runif(1,p,p+d)
  carr.mat.gamb[k,1]<-p+d
  carr.mat.gamb[k,2]<-t
}
carr.mat.gamb<-data.table(carr.mat.gamb)
colnames(carr.mat.gamb)<-c("duration","transmission")
# mean(carr.mat$transmission)
# quantile(carr.mat$transmission)
clear.mat[1,2]<-mean(carr.mat.gamb$duration,na.rm=T)
clear.mat[1,3]<-mean(carr.mat.gamb$transmission,na.rm=T)
clear.mat[1,4]<-sd(carr.mat.gamb$transmission,na.rm=T)
clear.mat[1,5:7]<-quantile(carr.mat.gamb$transmission,na.rm=T)[c(3,2,4)]



#########KILIFI KENYA, CHRISPIN CHAGUZA FRONTIERS ########
nboot=10000
carr.mat.kil<-matrix(nrow=nboot,ncol=2)
for(k in 1:nboot){
  durs<-rexp(1000,rate=0.032)
  durs[which(durs>maxCarriageDuration)]<-NA
  durs<-durs[!is.na(durs)]
  d<-sample(durs,1,replace = FALSE,prob=durs)
  p<-runif(1,1,5)####incubation period
  # t<-runif(1,0,p+d)
  t<-runif(1,p,p+d)
  carr.mat.kil[k,1]<-p+d
  carr.mat.kil[k,2]<-t
}
carr.mat.kil<-data.table(carr.mat.kil)
colnames(carr.mat.kil)<-c("duration","transmission")
# mean(carr.mat$transmission)
# quantile(carr.mat$transmission)

clear.mat[2,2]<-mean(carr.mat.kil$duration,na.rm=T)
clear.mat[2,3]<-mean(carr.mat.kil$transmission,na.rm=T)
clear.mat[2,4]<-sd(carr.mat.kil$transmission,na.rm=T)
clear.mat[2,5:7]<-quantile(carr.mat.kil$transmission,na.rm=T)[c(3,2,4)]


###PLOT#####
clear.mat<-data.frame(data.table(clear.mat))
clear.mat[,2:6]<-sapply(clear.mat[,2:6], as.numeric)
ggplot(clear.mat)+
  geom_point(aes(x=site,y=gentime_mean))+
  geom_errorbar(aes(ymin=gentime_mean-sd,ymax=gentime_mean+sd,x=site),width=0.2)+
  theme_bw()+
  ylab("Estimated Generation Time")+
  xlab("Site")

mean(clear.mat[,"gentime_mean"])
mean(clear.mat[,"sd"])

genTime<- mean(clear.mat[,"gentime_mean"])
varGen<-(genTime)^2
shape=genTime^2/varGen
scale=varGen/genTime

tots=250
colors=c("The Gambia"="grey","Kenya"="pink")
ggplot()+
  geom_histogram(data=carr.mat.kil[sample(tots,replace=F),], aes(x=transmission,y = ..count../sum(..count..),fill="grey"),binwidth=1,alpha=0.8)+
  geom_histogram(data=carr.mat.gamb[sample(tots,replace=F),], aes(x=transmission,y = ..count../sum(..count..),fill="pink"),binwidth=1,alpha=0.5,fill="pink")+
  geom_line(aes(x=1:tots,y=dgamma(1:tots,shape=shape,scale=scale)),linetype="dashed")+
  theme_classic()+
  xlab("Transmission Generation Interval (days)")+
  ylab("Density")+
  scale_x_continuous(breaks=c(0,100,200),limits=c(0,200),labels=c(0,100,200))+
  scale_y_continuous(breaks=c(0,.025,0.05),limits=c(0,0.05),labels=c(0,0.025,0.05))+
  # xlim(0,0.025)+
  theme(legend.title = element_blank(),legend.position = c(0.6,0.8),axis.text = element_text(size=15),axis.title = element_text(size=15))+
  scale_fill_manual(values = colors,breaks=c("The Gambia","Kenya"),limits=c("The Gambia","Kenya"))+
  guides(fill = guide_legend(override.aes = list(size=10,color="white")))

# ggsave("./CarriageDurations/transmissiongens.pdf",width=4,height=3.5)
  