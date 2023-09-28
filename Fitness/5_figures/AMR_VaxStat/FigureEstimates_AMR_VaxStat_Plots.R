## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(ggplot2)

setwd("./Fitness")
## Load fit

# fit = readRDS(file = "/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/RunModel/AMR_VaxStat/FixedTogether19JUN22/Output_vaxstatAMR_NVT.R_fit_all.rds")
fit = readRDS(file = "./4_run_model/AMR_VaxStat/output/Output_vaxstatAMR_NVT.R_fit_all.rds")


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
fitness_genotypes_pre_vacc_total[,1]<-Chains$fitness_genotypes_pre_vacc_amr[,1,1]+rep(0,nb_chains)
fitness_genotypes_pre_vacc_total[,2]<-rep(0,nb_chains)+rep(0,nb_chains)
fitness_genotypes_pre_vacc_total[,3]<-Chains$fitness_genotypes_pre_vacc_amr[,2,1]+Chains$fitness_genotypes_pre_vacc_vstat[,,1]
fitness_genotypes_pre_vacc_total[,4]<-rep(0,nb_chains)+Chains$fitness_genotypes_pre_vacc_vstat[,,1]

fitness_genotypes_post_vacc_total[,1]<-Chains$fitness_genotypes_post_vacc_amr[,1,1]+rep(0,nb_chains)
fitness_genotypes_post_vacc_total[,2]<-rep(0,nb_chains)+rep(0,nb_chains)
fitness_genotypes_post_vacc_total[,3]<-Chains$fitness_genotypes_post_vacc_amr[,2,1]+Chains$fitness_genotypes_post_vacc_vstat[,,1]
fitness_genotypes_post_vacc_total[,4]<-rep(0,nb_chains)+Chains$fitness_genotypes_post_vacc_vstat[,,1]

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
save(df_overall_fitness_AMRVaxStat,file="./5_figures/AMR_VaxStat/df_overall_fitness_AMRVaxStat.RData")
################################################################################
## Fitness, per serotype, per vaccine era (1 switch in 2010)
################################################################################
pdf(file = './5_figures/AMR_VaxStat/Figure_estimates_overall_switch2010.pdf', width = 9/2.54, height = 6/2.54)
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('Resistant','Susceptible'))
df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels =  c("nvt_s","nvt_r","vt_s","vt_r"))

################################################################################
ACV = ggplot(df_overall_fitness, aes(x = Time, y=Values, color = Genotype)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f2, lwd = 0.2, alpha=1, position = position_dodge2(width = 0.2), geom = "pointrange")+ 
  theme_classic()+
  facet_grid(.~Genotype)+
  # geom_hline(data = simple_model_R, aes(yintercept = Z))+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  # geom_vline(xintercept = -5, linetype = "longdash", color = 'grey')+
  scale_x_continuous(limits = c(-2,2), 
                     breaks = c(-1,1),
                     labels = c('2000-2009', '2010-2015'))+
  scale_y_continuous(trans = 'log', limits = c(0.6,1.6),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "Clades", values = c('Black', 'Black','Black', 'Black'))+
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 10, hjust=1, angle = 30),
        axis.text.y = element_text(size = 10,hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ACV
dev.off()
res = data.frame('group' = c('nvt_s','nvt_r','vt_s','vt_r'),
                 'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'nvt_s')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'nvt_r')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'vt_s')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'vt_r')]), digits = 2), collapse =  '-')),
                 'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'nvt_s')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'nvt_r')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'vt_s')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'vt_r')]), digits = 2), collapse =  '-')))
# res_mean = data.frame('group' = c('nvt_s','nvt_r','vt_s','vt_r'),
#                       'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Resistant')])[1], digits = 2), collapse =  '-'),
#                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Susceptible')])[1], digits = 2), collapse =  '-')),
#                       'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Resistant')])[1], digits = 2), collapse =  '-'),
#                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Susceptible')])[1], digits = 2), collapse =  '-')))
# res_cimin = data.frame('group' = c('nvt_s','nvt_r','vt_s','vt_r'),
#                       'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Resistant')])[2], digits = 2), collapse =  '-'),
#                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Susceptible')])[2], digits = 2), collapse =  '-')),
#                       'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Resistant')])[2], digits = 2), collapse =  '-'),
#                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Susceptible')])[2], digits = 2), collapse =  '-')))
# res_cimax = data.frame('group' = c('nvt_s','nvt_r','vt_s','vt_r'),
#                        'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Resistant')])[3], digits = 2), collapse =  '-'),
#                                                paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Susceptible')])[3], digits = 2), collapse =  '-')),
#                        'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Resistant')])[3], digits = 2), collapse =  '-'),
#                                                paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Susceptible')])[3], digits = 2), collapse =  '-')))
saveRDS(res, './4_run_model/AMR_VaxStat/output/Results_output_per_provice_pooled_by_vaccine_swicth2009_plus1.rds')

################################################################################

################################################################################
# ## Fitness, per serotype, per vaccine era (2 switches in 2009 and 2011)
# ################################################################################
# # pdf(file = 'Figure_estimates_switch2009_2011.pdf', width = 12/2.54, height = 7/2.54)
# # windows(width = 15/2.54, height = 8/2.54)
# ################################################################################
# df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = 1,
#                                                           'Time' = rep(-3, each = nb_chains),
#                                                           'Genotype' = rep(Genotype[1], nb_chains*length_vect_fitness_post_vacc)))
# df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = 1,
#                                                           'Time' = rep(-3, each = nb_chains),
#                                                           'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('Resistant' ,'Susceptible'))
# ACV = ggplot(df_overall_fitness, aes(x = Time, y=Values, color = Genotype)) + 
#   geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
#   stat_summary(fun.data = f2, lwd = 0.2, alpha=1, position = position_dodge2(width = 0.2), geom = "pointrange")+ 
#   theme_classic()+
#   facet_grid(.~Genotype)+
#   # facet_grid(.~factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13')))+
#   # geom_hline(data = simple_model_R, aes(yintercept = Z))+
#   # geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
#   # geom_vline(xintercept = -2, linetype = "longdash", color = 'grey')+
#   # geom_vline(xintercept = -5, linetype = "longdash", color = 'grey')+
#   scale_x_continuous(limits = c(-4,2), 
#                      breaks = c(-3, -1 ,1),
#                      labels = c('2000-2008', '2009-2010', '2011-2015'))+
#   scale_y_continuous(trans = 'log', limits = c(0.6,1.6),
#                      breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
#                      labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
#   labs(title = "", x = '', y = 'Fitness')+
#   scale_color_manual(name = "Clades", values = c("red", "forestgreen"))+                                                  
#   theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10),
#         axis.title.x = element_text(size=10),
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size = 10, hjust=1, angle = 30),
#         axis.text.y = element_text(size = 10, hjust = 0.5),
#         strip.text.x = element_text(size = 10, colour = "black", angle = 0),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")
# 
# ACV
# # dev.off()
# res = data.frame('PENR' = c('Resistant','Susceptible'),
#                  'Fitness_2009-2010' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Resistant')]), digits = 2), collapse =  '-'),
#                                          paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Susceptible')]), digits = 2), collapse =  '-')),
#                  'Fitness_2011-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Resistant')]), digits = 2), collapse =  '-'),
#                                          paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Susceptible')]), digits = 2), collapse =  '-')))

################################################################################




