####Figure 4 D-G Fitness Estimates pre- post- and for each group
library(ggplot2)
library(rstan)

df_overall_fitness<-readRDS("./4_run_model/NVT_PCV7_PCV13/output/Data_plot_per_VT_22012024.rds")
## Load fit
# fit = readRDS(file = './4_run_model/NVT_PCV7_PCV13/output/Output_per_provice_NVT_PCV7_PCV13_swicth2009_plus0_fit_all.rds')
fit = readRDS(file = './4_run_model/NVT_PCV7_PCV13/output/Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011_1pervax_plus0_fit_all.rds')

## Chains
Chains=rstan::extract(fit$fit)
#######FUNCTIONS
mean.and.ci <-function(v){ 
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}
f2 <- function(x) {
  r <- quantile(x, probs = c(0.025,0.5,0.975))
  names(r) <- c("ymin","y","ymax")
  r
}


df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
####Plot pre and post
p1<-ggplot(df_overall_fitness,aes(x=Time,y=Values,color=Genotype))+
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f2, lwd = 0.2, alpha=1, position = position_dodge2(width = 0.2), geom = "pointrange")+
  scale_x_continuous(limits = c(-2,2), 
                     breaks = c(-1,1),
                     labels = c('2000-2009', '2010-2015'))+
  scale_y_continuous(trans = 'log', limits = c(0.6,1.6),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+     
  theme_classic()+
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 10, hjust=1, angle = 30),
        axis.text.y = element_text(size = 10,hjust = 0.5)  )

################################################################################
## Compute VT fitness estimate - PAPER ESTIMATES
################################################################################
# fitness_VT_pre = rep(0, length(nb_chains))
# fitness_VT_post = rep(0, length(nb_chains))
# mean_freq_PCV7_pre = NULL
# mean_freq_PCV7_post = NULL
# mean_freq_PCV13_pre = NULL
# mean_freq_PCV13_post = NULL
# for(k in 1:nb_countries){
#   mean_freq = Chains$pred_absolute_freq[,k,1,vacc_WCV_k:(vacc_ACV_k-1)]
#   mean_freq = apply(mean_freq, MARGIN = 1, mean)
#   mean_freq_PCV7_pre = c(mean_freq_PCV7_pre, mean_freq)
#   
#   mean_freq = Chains$pred_absolute_freq[,k,1,vacc_ACV_k:nb_years]
#   mean_freq = apply(mean_freq, MARGIN = 1, mean)
#   mean_freq_PCV7_post = c(mean_freq_PCV7_post, mean_freq)
#   
#   mean_freq = Chains$pred_absolute_freq[,k,2,vacc_WCV_k:(vacc_ACV_k-1)]
#   mean_freq = apply(mean_freq, MARGIN = 1, mean)
#   mean_freq_PCV13_pre = c(mean_freq_PCV13_pre, mean_freq)
#   
#   mean_freq = Chains$pred_absolute_freq[,k,2,vacc_ACV_k:nb_years]
#   mean_freq = apply(mean_freq, MARGIN = 1, mean)
#   mean_freq_PCV13_post = c(mean_freq_PCV13_post, mean_freq)
# }
# fitness_VT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre) +
#   df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre)
# 
# fitness_VT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_post/(mean_freq_PCV7_post+mean_freq_PCV13_post) +
#   df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_post/(mean_freq_PCV7_post+mean_freq_PCV13_post)




fitness_NVT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')]
fitness_NVT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')]
fitness_PCV7_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')]
fitness_PCV7_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')]
fitness_PCV13_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')]
fitness_PCV13_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')]

mean.and.ci(fitness_PCV7_post/fitness_PCV7_pre)
mean.and.ci(fitness_PCV13_post/fitness_PCV13_pre)
mean.and.ci(fitness_NVT_post/fitness_NVT_pre)

######Pre and post compared to an NVT reference####
prePostNVTRef.df<-matrix(nrow=6,ncol=5)
prePostNVTRef.df[1,1:3]<-mean.and.ci(fitness_NVT_pre/fitness_NVT_pre)
prePostNVTRef.df[2,1:3]<-mean.and.ci(fitness_PCV7_pre/fitness_NVT_pre)
prePostNVTRef.df[3,1:3]<-mean.and.ci(fitness_PCV13_pre/fitness_NVT_pre)
prePostNVTRef.df[1:3,4]<-"Pre"
prePostNVTRef.df[1,5]<-"NVT"
prePostNVTRef.df[2,5]<-"PCV7"
prePostNVTRef.df[3,5]<-"PCV13"
prePostNVTRef.df[4,1:3]<-mean.and.ci(fitness_NVT_post/fitness_NVT_post)
prePostNVTRef.df[5,1:3]<-mean.and.ci(fitness_PCV7_post/fitness_NVT_post)
prePostNVTRef.df[6,1:3]<-mean.and.ci(fitness_PCV13_post/fitness_NVT_post)
prePostNVTRef.df[4:6,4]<-"Post"
prePostNVTRef.df[4,5]<-"NVT"
prePostNVTRef.df[5,5]<-"PCV7"
prePostNVTRef.df[6,5]<-"PCV13"
prePostNVTRef.df<-data.frame(prePostNVTRef.df)

cols.num<-c("X1","X2","X3")
prePostNVTRef.df[cols.num]<-sapply(prePostNVTRef.df[cols.num],as.character)
prePostNVTRef.df[cols.num]<-sapply(prePostNVTRef.df[cols.num],as.numeric)
colnames(prePostNVTRef.df)<-c("mean","lowerCI","upperCI","Time","Group")
prePostNVTRef.df$Time<-factor(prePostNVTRef.df$Time,levels=c("Pre","Post"),labels=c("Pre-PCV","Post-PCV"))
prePostNVTRef.df$Group<-factor(prePostNVTRef.df$Group,levels=c("NVT","PCV7","PCV13"))


######Post compared to an pre as reference####
prePost.df<-matrix(nrow=6,ncol=5)
prePost.df[1,1:3]<-mean.and.ci(fitness_NVT_pre/fitness_NVT_pre)
prePost.df[2,1:3]<-mean.and.ci(fitness_PCV7_pre/fitness_PCV7_pre)
prePost.df[3,1:3]<-mean.and.ci(fitness_PCV13_pre/fitness_PCV13_pre)
prePost.df[1:3,4]<-"Pre"
prePost.df[1,5]<-"NVT"
prePost.df[2,5]<-"PCV7"
prePost.df[3,5]<-"PCV13"
prePost.df[4,1:3]<-mean.and.ci(fitness_NVT_post/fitness_NVT_pre)
prePost.df[5,1:3]<-mean.and.ci(fitness_PCV7_post/fitness_PCV7_pre)
prePost.df[6,1:3]<-mean.and.ci(fitness_PCV13_post/fitness_PCV13_pre)
prePost.df[4:6,4]<-"Post"
prePost.df[4,5]<-"NVT"
prePost.df[5,5]<-"PCV7"
prePost.df[6,5]<-"PCV13"
prePost.df<-data.frame(prePost.df)
cols.num<-c("X1","X2","X3")
prePost.df[cols.num]<-sapply(prePost.df[cols.num],as.character)
prePost.df[cols.num]<-sapply(prePost.df[cols.num],as.numeric)
colnames(prePost.df)<-c("mean","lowerCI","upperCI","Time","Group")
prePost.df$Time<-factor(prePost.df$Time,levels=c("Pre","Post"),labels=c("Pre-PCV","Post-PCV"))
prePost.df$Group<-factor(prePost.df$Group,levels=c("NVT","PCV7","PCV13"))
# 
prePost.df.estimates<-prePost.df
save(prePost.df.estimates,file="./5_figures/AMR_VaxStat/prePost.df.estimates.RData")
prePostNVTRef.df.estimates<-prePostNVTRef.df
save(prePostNVTRef.df.estimates,file="./5_figures/AMR_VaxStat/prePostNVTRef.df.estimates.RData")

p2<-ggplot(prePostNVTRef.df.estimates,aes(x=Group,y=mean,group=Time,color=Group))+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=5,shape=18)+
  theme_classic()+
  geom_errorbar(data=subset(prePostNVTRef.df.estimates,prePostNVTRef.df.estimates$Group!="NVT"),aes(x=Group,ymin=lowerCI,ymax=upperCI,group=Time,color=Group),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.5,1.6),
                     breaks = c(0.1, 0.25,0.5, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text=element_text(size=20),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~Time)

p3<-ggplot(prePost.df.estimates,aes(x=Time,y=mean,group=Group,color=Group))+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=3,shape=18)+
  theme_classic()+
  geom_errorbar(data=subset(prePost.df.estimates,prePost.df.estimates$Time!="Pre-PCV"), aes(x=Time,ymin=lowerCI,ymax=upperCI,group=Time,color=Group),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.5,1.6),
                     breaks = c(0.1, 0.25,0.5,  1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text.x=element_text(size=20,angle=45,vjust=0.5),axis.text=element_text(size=20),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~Group)

library(patchwork)
p2/p3
p3

load('./5_figures/AMR_VaxStat/df_overall_fitness_AMRVaxStat.RData')
df_overall_fitness_AMRVaxStat
fitness_nvtR_pre = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '-1' & df_overall_fitness_AMRVaxStat$Genotype == 'nvt_r')]
fitness_nvtR_post = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '1' & df_overall_fitness_AMRVaxStat$Genotype == 'nvt_r')]
fitness_nvtS_pre = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '-1' & df_overall_fitness_AMRVaxStat$Genotype == 'nvt_s')]
fitness_nvtS_post = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '1' & df_overall_fitness_AMRVaxStat$Genotype == 'nvt_s')]
fitness_vtR_pre = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '-1' & df_overall_fitness_AMRVaxStat$Genotype == 'vt_r')]
fitness_vtR_post = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '1' & df_overall_fitness_AMRVaxStat$Genotype == 'vt_r')]
fitness_vtS_pre = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '-1' & df_overall_fitness_AMRVaxStat$Genotype == 'vt_s')]
fitness_vtS_post = df_overall_fitness_AMRVaxStat$Values[which(df_overall_fitness_AMRVaxStat$Time == '1' & df_overall_fitness_AMRVaxStat$Genotype == 'vt_s')]

mean.and.ci(fitness_nvtR_post/fitness_nvtR_pre)
mean.and.ci(fitness_nvtR_post/fitness_nvtS_post)
mean.and.ci(fitness_vtR_post/fitness_vtR_pre)
mean.and.ci(fitness_vtR_post/fitness_vtS_post)

mean.and.ci(fitness_nvtR_post+fitness_vtR_post/2)



#### numbers in the paper
mean.and.ci(fitness_NVT_pre/c(fitness_PCV7_pre,fitness_PCV13_pre))
mean.and.ci(c(fitness_PCV7_pre,fitness_PCV13_pre)/fitness_NVT_pre)
mean.and.ci(c(fitness_PCV7_post,fitness_PCV13_post)/fitness_NVT_post)
mean.and.ci(c(fitness_PCV7_post)/fitness_NVT_post)
mean.and.ci(c(fitness_PCV13_post)/fitness_NVT_post)


mean.and.ci(fitness_NVT_pre)
mean.and.ci(fitness_NVT_post)

mean.and.ci(fitness_PCV7_pre)
mean.and.ci(fitness_PCV7_post)

mean.and.ci(fitness_PCV13_pre)
mean.and.ci(fitness_PCV13_post)



