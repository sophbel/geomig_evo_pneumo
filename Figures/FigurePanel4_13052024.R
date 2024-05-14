#############
######FIGURE 4##########
################
library(ggplot2)
library(data.table)
##### Figure 4A-C#####------------------------------------------------------------
load("./Figures/Figure4/fit_table_nvtpcv7pcv13_2switch.RData")
mat.all<-mat_nvtvt
mat.all$type_f<-factor(mat.all$type,levels=c("NVT","PCV7","PCV13"))
AC4<-ggplot()+
  geom_point(data=mat.all,aes(x=years,y=data,group=type_f,color=type_f))+
  geom_line(data=mat.all,aes(x=years,y=fit,group=type_f,color=type_f),alpha=0.8)+
  theme_classic()+
  geom_errorbar(data=mat.all,aes(x=years,ymin=datalower,ymax=dataupper,group=type_f,color=type_f))+
  geom_ribbon(data=mat.all,aes(x=years,ymin=fitlower,ymax=fitupper,group=type_f,fill=type_f),alpha=0.5)+
  scale_color_manual(values = c("NVT"="maroon","PCV7"="darkblue","PCV13"="darkgreen"))+
  scale_fill_manual(values = c("NVT"="maroon","PCV7"="darkblue","PCV13"="darkgreen"))+
  theme(axis.text.x=element_text(size=20,angle=45,vjust=0.6),axis.text.y=element_text(size=20),axis.title = element_text(size=20),plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "none")+
  ylim(0,1)+
  xlab("Years")+
  ylab("Proportion")+
  facet_wrap(type_f~.)+
  theme(panel.spacing = unit(2, "lines"),strip.text.x = element_text(size=18))
AC4


#####Figure 4D-E######------------------------------------------------------------
### mat ests relative fitness
mat_ests<-fread("./Figures/Figure4/Estimates_Fig4DE_relative_fitness_estimates_prevax_NVT_PCV7_PCV13_22012026.csv",header = TRUE)
mat_ests$Time[which(mat_ests$Time=="Pre.PCV")]<-c("Relative Fitness Pre-PCV")
mat_ests$Time[which(mat_ests$Time=="Post.PCV")]<-c("Relative Fitness Post-PCV")
mat_ests$Time<-as.factor(mat_ests$Time)
mat_ests$Time<-factor(mat_ests$Time,level=c("Relative Fitness Pre-PCV","Relative Fitness Post-PCV"))
mat_ests$Group<-factor(mat_ests$Group,level=c("NVT","PCV7","PCV13"))
### pre pcv as reference new
### 
mat_ests_preref<-fread("./Figures/Figure4/Estimates_Fig4F_relative_fitness_estimates_post_vs_prevax_NVT_PCV7_PCV13_22012026.csv",header = TRUE)
mat_ests_preref$Group<-as.factor(mat_ests_preref$Group)
mat_ests_preref$Group<-factor(mat_ests_preref$Group,level=c("NVT","PCV7","PCV13"))
mat_ests_preref$Time<-factor(mat_ests_preref$Time,level=c("Pre-PCV","Post-PCV"))
### 
D4<-ggplot(mat_ests,aes(x=Group,y=mean,group=Time,color=Group))+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=5,shape=18)+
  theme_classic()+
  geom_errorbar(data=subset(mat_ests,mat_ests$Group!="NVT"),aes(x=Group,ymin=lowerCI,ymax=upperCI,group=Time,color=Group),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.49,1.6),
                     breaks = c(0.1, 0.25,0.5, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text=element_text(size=20),axis.text.x=element_text(vjust=0.6),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~Time)

E4<-ggplot(mat_ests_preref,aes(x=Time,y=mean,group=Group,color=Group))+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=5,shape=18)+
  theme_classic()+
  geom_errorbar(data=subset(mat_ests_preref,mat_ests_preref$Time!="Pre-PCV"), aes(x=Time,ymin=lowerCI,ymax=upperCI,group=Time,color=Group),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.5,1.6),
                     breaks = c(0.1, 0.25,0.5,  1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text.x=element_text(size=20,angle=45,vjust=0.5),axis.text=element_text(size=20),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~Group)

library(patchwork)
D4+E4


##### Figure 4F#####------------------------------------------------------------
# Combining the individual resistance fits into the overall plot ggplot
load("./Figures/Figure4/df.projection.RData")###projectR_vaxStat_newModel
df.projection<-subset(df.projection,df.projection$years<2016)
df.projection[which(df.projection$years==2015),]<-NA
load("./Figures/Figure4/df.vaxstatamr.newmodel.RData")####AMR_VTNVTtogether_justAMR

df.vaxstatamr<-data.frame(df.vaxstatamr)
df.vaxstatamr[,c(1:7)]<-apply(df.vaxstatamr[,c(1:7)],2,function(x) as.numeric(x))
df.vaxstatamr$years<-df.vaxstatamr$years+1999
df.vaxstatamr.R<-subset(df.vaxstatamr,df.vaxstatamr$amr=="Resistant")
#

F4=ggplot()+
  geom_point(data=df.projection,aes(x=years,y=data),color="black")+
  geom_point( data=df.vaxstatamr.R,aes(x=years,y=data,color=vaxstat),show.legend = T)+
  geom_errorbar(data=df.vaxstatamr.R,aes(x=years,ymin=di_lower,ymax=di_upper,color=vaxstat),width=0.2,show.legend = F)+
  geom_ribbon(data=df.vaxstatamr.R,aes(x=years,ymin=fit_lower,ymax=fit_upper,fill=vaxstat),color=NA,alpha=0.6)+
  geom_ribbon(data=df.projection,aes(x=years,ymin=fit_lower,ymax=fit_upper),fill="black",color=NA,alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = c("maroon","darkblue"))+
  scale_fill_manual(values = c("maroon","darkblue"))+
  scale_x_discrete(breaks=c(2000,2005,2010,2015),limits=c(2000,2005,2010,2015),labels=c(2000,2005,2010,2015))+
  ylim(0,1)+
  theme(legend.title = element_blank(),legend.position = c(10,0.9),axis.text=element_text(size=20),axis.text.x = element_text(angle=45,vjust=0.6),
        axis.title = element_text(size=20),
        legend.text = element_text(size=15),
        plot.margin = margin(0.5, 0.5,0.5, 0.5, "cm"))+
  ylab("Proportion Resistant")+
  xlab("Years")
F4

#####Figure 4G#####------------------------------------------------------------
fit_R_main<-read.csv(file="./Figures/Figure4/fitRmain_4G.csv")
fit_R_main$time<-factor(fit_R_main$time,levels=c("Pre-PCV","Post-PCV"))
G4<-ggplot(fit_R_main,aes(x=type,y=mean,color=type,group=time))+
  geom_point()+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=5,shape=18)+
  theme_classic()+
  geom_errorbar(data=fit_R_main, aes(x=type,ymin=lowerCI,ymax=upperCI,color=type,group=time),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.8,1.6), breaks = c(0.8,  1.0, 1.5, 2, 3, 10), labels = c(0.8,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative Fitness of AMR')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue'))+
  theme(axis.text.x=element_text(size=20,angle=45,vjust=0.5),axis.text=element_text(size=20),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~time)
G4


library(patchwork)
AC4/ (D4+E4)/(F4+G4)


