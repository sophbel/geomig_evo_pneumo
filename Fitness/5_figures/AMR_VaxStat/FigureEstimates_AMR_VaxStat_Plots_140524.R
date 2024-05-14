## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(ggplot2)

# setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/RunModel/AMR_VaxStat/FixedTogether19JUN22/")
## Load fit

fit = readRDS(file = './4_run_model/AMR_VaxStat/output/Output_vaxstatAMR_NVT.R_fit_all.rds')


## Extract chains
Chains=rstan::extract(fit$fit)
################################################################################
## Functions
################################################################################
mean.and.ci <-function(v){ 
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}
# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
f2 <- function(x) {
  r <- quantile(x, probs = c(0.025,0.5,0.975))
  names(r) <- c("ymin","y","ymax")
  r
}
################################################################################

################################################################################
## Compute data to plot
################################################################################

################################################################################
## Growth rate (fitness) weighted by proportions
################################################################################
vaccination = fit$data$vaccine_introduction
nb_chains = length(Chains$lp__)
nb_countries = fit$data$nb_countries
nb_years = fit$data$nb_years
length_vect_fitness_post_vacc = 1#fit$data$number_R_post_vacc
length_vect_fitness_pre_vacc = 1#fit$data$number_R_pre_vacc
nb_genotypes_total = fit$data$nb_genotypes_total

ref_genotype = 2 ###NVT_R - 2nd row


ref_genotype_header = "nvt_r"
Genotype = c("nvt_s","vt_s","vt_r")

####Total Fitnesses
# vt.type <- c("NVT","VT") ref_clade_vstat = 1 ## Ref clade: NVT
# amr.type <- c(0,1) ref_clade_amr = 2 ## Ref clade: Resistant

fitness_genotypes_pre_vacc_total<-matrix(nrow=nb_chains,ncol=nb_genotypes_total)
fitness_genotypes_post_vacc_total<-matrix(nrow=nb_chains,ncol=nb_genotypes_total)

###Colnames("nvt_s","nvt_r","vt_s","vt_r")
fitness_genotypes_pre_vacc_total[,1]<-Chains$fitness_genotypes_pre_vacc_amr[,1,1]+rep(0,nb_chains)# non vaccine type susc
fitness_genotypes_pre_vacc_total[,2]<-rep(0,nb_chains)+rep(0,nb_chains)
fitness_genotypes_pre_vacc_total[,3]<-Chains$fitness_genotypes_pre_vacc_amr[,2,1]+Chains$fitness_genotypes_pre_vacc_vstat[,,1]
fitness_genotypes_pre_vacc_total[,4]<-rep(0,nb_chains)+Chains$fitness_genotypes_pre_vacc_vstat[,,1]

fitness_genotypes_post_vacc_total[,1]<-Chains$fitness_genotypes_post_vacc_amr[,1,1]+rep(0,nb_chains)
fitness_genotypes_post_vacc_total[,2]<-rep(0,nb_chains)+rep(0,nb_chains)
fitness_genotypes_post_vacc_total[,3]<-Chains$fitness_genotypes_post_vacc_amr[,2,1]+Chains$fitness_genotypes_post_vacc_vstat[,,1]
fitness_genotypes_post_vacc_total[,4]<-rep(0,nb_chains)+Chains$fitness_genotypes_post_vacc_vstat[,,1]

## ---------------------------------------------------------------
##### calculate f0 proportions for each in 2009. 
  ###Colnames("nvt_s","nvt_r","vt_s","vt_r")
  
###NVT_R
mean.and.ci(apply(Chains$pred_absolute_freq_total[,,2,10],1,mean))
###NVT_S
mean.and.ci(apply(Chains$pred_absolute_freq_total[,,1,10],1,mean))




##-------------------------------------------------------------------------------------------------------------------------------------

# Chains$fitness_genotypes_pre_vacc_amr ###S_NVT and S_VT
# Chains$fitness_genotypes_pre_vacc_vstat ### VT Column
# 
# Chains$fitness_genotypes_post_vacc_amr ###S_NVT and S_VT
# Chains$fitness_genotypes_post_vacc_vstat ###VT Column


# ref_genotype_header = '1'
# Genotype = rownames(clade_sero_list[[1]])[-1]
df_overall_fitness = data.frame('Values' = NA,
                                'Time' = NA,
                                'Genotype' = NA)
vaccinationACV = fit$data$yearIntroduction
vaccinationWCV = fit$data$yearF0

Genotype<-c("nvt_s","nvt_r","vt_s","vt_r")
for(k in 1:nb_countries){
  for(g in 1:(nb_genotypes_total)){
    
    ## Compute reference fitness first:
    ## pre vacc
    mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes_total, nrow = nb_chains)
    vacc_ACV_k = vaccinationACV[k]
    vacc_WCV_k = vaccinationWCV[k]
    for(gref in 1:(nb_genotypes_total)){
      mat_tmp[(1:nb_chains),gref] = as.vector(fitness_genotypes_pre_vacc_total[,g])-as.vector(fitness_genotypes_pre_vacc_total[,gref])
      mean_freq = Chains$pred_absolute_freq_total[,k,gref,vacc_WCV_k:(vacc_ACV_k-1)]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains),gref] = mat_tmp[(1:nb_chains),gref] * mean_freq
    }
    
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                              'Time' = rep(-1, each = nb_chains),
                                                              'Genotype' = rep(Genotype[g], nb_chains)))
    
    ## post vacc
    mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes_total, nrow = nb_chains)
    for(gref in 1:(nb_genotypes_total)){
      mat_tmp[(1:nb_chains),gref] = as.vector(fitness_genotypes_post_vacc_total[,g])-as.vector(fitness_genotypes_post_vacc_total[,gref])
      mean_freq = Chains$pred_absolute_freq_total[,k,gref,vacc_ACV_k:nb_years]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains),gref] = mat_tmp[(1:nb_chains),gref] * mean_freq
    }
    
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                              'Time' = rep(1, each = nb_chains),
                                                              'Genotype' = rep(Genotype[g], nb_chains)))
  }
}

df_overall_fitness = df_overall_fitness[-1,] ## Remove the NA from the beginning


################################################################################
## Set directory
################################################################################
df_overall_fitness_AMRVaxStat<-df_overall_fitness
# save(df_overall_fitness_AMRVaxStat,file="/Users/sophiebelman/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/Fitness/testing_AMRVaxStat/df_overall_fitness_AMRVaxStat.RData")
################################################################################
## Fitness, per serotype, per vaccine era (1 switch in 2010)
################################################################################
# pdf(file = '/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/RunModel/AMR_VaxStat/Fits_Estimates/Figure_estimates_overall_switch2010.pdf', width = 9/2.54, height = 6/2.54)
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('Resistant','Susceptible'))
df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels =  c("nvt_s","nvt_r","vt_s","vt_r"))




### compute to plot

# levels =  c("nvt_s","nvt_r","vt_s","vt_r")
fitness_NVT_R_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'nvt_r')]
fitness_NVT_R_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'nvt_r')]

fitness_NVT_S_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'nvt_s')]
fitness_NVT_S_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'nvt_s')]

fitness_VT_R_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'vt_r')]
fitness_VT_R_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'vt_r')]

fitness_VT_S_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'vt_s')]
fitness_VT_S_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'vt_s')]

mean.and.ci(fitness_NVT_R_pre)

mean.and.ci(fitness_NVT_S_pre)
mean.and.ci(fitness_NVT_S_post)
mean.and.ci(fitness_NVT_R_post)
mean.and.ci(fitness_VT_R_post)
mean.and.ci(fitness_VT_S_post)

#### relative fitness post to pre PCV
groups<-c("NVT_S","NVT_R","VT_S","VT_R")
fit_vaxstat<-matrix(1,nrow=8,ncol = 5)
fit_vaxstat[1,2:4]<-mean.and.ci(fitness_NVT_S_post/fitness_NVT_S_pre)
fit_vaxstat[2,2:4]<-mean.and.ci(fitness_NVT_R_post/fitness_NVT_R_pre)
fit_vaxstat[3,2:4]<-mean.and.ci(fitness_VT_S_post/fitness_VT_S_pre)
fit_vaxstat[4,2:4]<-mean.and.ci(fitness_VT_R_post/fitness_VT_R_pre)
fit_vaxstat[,1]<-rep(groups,2)
fit_vaxstat[,5]<-c(rep("Post-PCV",4),rep("Pre-PCV",4))
fit_vaxstat<-data.frame(fit_vaxstat)
fit_vaxstat[,2:4]<-(sapply(fit_vaxstat[,2:4],as.numeric))
colnames(fit_vaxstat)<-c("group","mean","lowerCI","upperCI","time")
################################################################################

fit_vaxstat$time<-factor(fit_vaxstat$time,levels=c("Pre-PCV","Post-PCV"))
fit_vaxstat$group<-factor(fit_vaxstat$group,levels=c("NVT_R","NVT_S","VT_R","VT_S"))
(p<-ggplot(fit_vaxstat,aes(x=time,y=mean,group=group,color=group))+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  geom_point(size=5,shape=18)+
  theme_classic()+
  geom_errorbar(data=fit_vaxstat, aes(x=time,ymin=lowerCI,ymax=upperCI,group=time,color=group),width=0.3)+
  scale_y_continuous(trans = 'log', limits = c(0.8,1.6),
                     breaks = c(0.8,  1.0, 1.5, 2, 3, 10), 
                     labels = c('0.8',   1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkgreen', 'gold','darkblue'))+
  theme(axis.text.x=element_text(size=20,angle=45,vjust=0.5),axis.text=element_text(size=20),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
  facet_grid(.~group))
# ggsave(p,file="/Users/sophiebelman/Documents/Migration/SA_Migration_110422/Submission/Submission2/Nature_Revisions/Revisions3_260324/Estimates_AMRVaxStat/Fit_Estimates_S22_030424.pdf",width=6.3,height=4)





#### absolute fitness estimates

groups<-c("NVT_S","NVT_R","VT_S","VT_R")
fit_vaxstat_abs<-matrix(1,nrow=8,ncol = 5)
fit_vaxstat_abs[1,2:4]<-mean.and.ci(fitness_NVT_S_pre)
fit_vaxstat_abs[2,2:4]<-mean.and.ci(fitness_NVT_R_pre)
fit_vaxstat_abs[3,2:4]<-mean.and.ci(fitness_VT_S_pre)
fit_vaxstat_abs[4,2:4]<-mean.and.ci(fitness_VT_R_pre)
fit_vaxstat_abs[5,2:4]<-mean.and.ci(fitness_NVT_S_post)
fit_vaxstat_abs[6,2:4]<-mean.and.ci(fitness_NVT_R_post)
fit_vaxstat_abs[7,2:4]<-mean.and.ci(fitness_VT_S_post)
fit_vaxstat_abs[8,2:4]<-mean.and.ci(fitness_VT_R_post)
fit_vaxstat_abs[,1]<-rep(groups,2)
fit_vaxstat_abs[,5]<-c(rep("Pre-PCV",4),rep("Post-PCV",4))
fit_vaxstat_abs<-data.frame(fit_vaxstat_abs)
fit_vaxstat_abs[,2:4]<-(sapply(fit_vaxstat_abs[,2:4],as.numeric))
colnames(fit_vaxstat_abs)<-c("group","mean","lowerCI","upperCI","time")
################################################################################

fit_vaxstat_abs$time<-factor(fit_vaxstat_abs$time,levels=c("Pre-PCV","Post-PCV"))
fit_vaxstat_abs$group<-factor(fit_vaxstat_abs$group,levels=c("NVT_R","NVT_S","VT_R","VT_S"))
(p_abs<-ggplot(fit_vaxstat_abs,aes(x=time,y=mean,group=group,color=group))+
    geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
    geom_point(size=5,shape=18)+
    theme_classic()+
    geom_errorbar(data=fit_vaxstat_abs, aes(x=time,ymin=lowerCI,ymax=upperCI,group=time,color=group),width=0.3)+
    scale_y_continuous(trans = 'log', limits = c(0.8,1.6),
                       breaks = c(0.8,  1.0, 1.5, 2, 3, 10), 
                       labels = c('0.8',   1.0, 1.5, 2, 3, 10))+
    labs(title = "", x = '', y = 'Relative Fitness')+
    scale_color_manual(name = "", values = c( 'maroon','darkgreen', 'gold','darkblue'))+
    theme(axis.text.x=element_text(size=20,angle=45,vjust=0.5),axis.text=element_text(size=20),axis.title = element_text(size=20),
          strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")+
    facet_grid(.~group))
# ggsave(p_abs,file="/Users/sophiebelman/Documents/Migration/SA_Migration_110422/Submission/Submission2/Nature_Revisions/Revisions3_260324/Estimates_AMRVaxStat/Fit_Estimates_S22_030424_absvals.pdf",width=6.3,height=4)

#### AMR fitness over susceptible
fit_R_main<-matrix(nrow=4,ncol=6)
fit_R_main[1,1:3]<-mean.and.ci(fitness_NVT_R_pre/fitness_NVT_S_pre)
fit_R_main[2,1:3]<-mean.and.ci(fitness_NVT_R_post/fitness_NVT_S_post)
fit_R_main[3,1:3]<-mean.and.ci(fitness_VT_R_pre/fitness_VT_S_pre)
fit_R_main[4,1:3]<-mean.and.ci(fitness_VT_R_post/fitness_VT_S_post)
fit_R_main[,4]<-c("NVT_pre","NVT_post","VT_pre","VT_post")
fit_R_main[,5]<-rep(c("Pre-PCV","Post-PCV"),2)
fit_R_main[,6]<-c(rep("NVT",2),rep("VT",2))


fit_R_main<-data.frame(fit_R_main)
fit_R_main[,1:3]<-(sapply(fit_R_main[,1:3],as.numeric))
colnames(fit_R_main)<-c("mean","lowerCI","upperCI","group","time","type")
fit_R_main$time<-factor(fit_R_main$time,levels=c("Pre-PCV","Post-PCV"))
fit_R_main$type<-factor(fit_R_main$type,levels=c("NVT","VT"))
tmp<-subset(fit_R_main,fit_R_main$time=="Post-PCV")

# write.csv(fit_R_main,file="/Users/sophiebelman/Documents/Migration/SA_Migration_110422/Submission/Submission2/Nature_Revisions/Final_Revisions_090524/Figures/Figure4/fitRmain_4G.csv",quote=FALSE,row.names = FALSE)
fit_R_main<-read.csv(file="./Fitness/5_figures/AMR_VaxStat/fitRmain_4G.csv")
fit_R_main$time<-factor(fit_R_main$time,levels=c("Pre-PCV","Post-PCV"))

pamr<-ggplot(fit_R_main,aes(x=type,y=mean,color=type,group=time))+
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

# ggsave(pamr,file="/Users/sophiebelman/Documents/Migration/SA_Migration_110422/Submission/Submission2/Nature_Revisions/Revisions3_260324/Estimates_AMRVaxStat/Fit_Estimates_AMR_Fig4_030424.pdf",width=6.3,height=4)



