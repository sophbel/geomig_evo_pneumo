## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)

# setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/')

## Load fit
fit = readRDS(file = './4_run_model/NVT_PCV7_PCV13/output/4_1_per_province_NVT_PCV7_PCV13_testswicth/disease/Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011_1pervax_disease_plus0_fit_all.rds')

## Chains
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
vaccination = fit$data$booster_introduction
nb_chains = length(Chains$lp__)
nb_countries = fit$data$nb_countries
nb_years = fit$data$nb_years
length_vect_fitness_post_vacc = fit$data$number_R_post_vacc
length_vect_fitness_pre_vacc = fit$data$number_R_pre_vacc
nb_genotypes = fit$data$nb_genotypes
ref_genotype = 1
ref_genotype_header = 'NVT'
Genotype = c('PCV7', 'PCV13')

# ref_genotype_header = '1'
# Genotype = rownames(clade_sero_list[[1]])[-1]

df_overall_fitness = data.frame('Values' = NA,
                                'Time' = NA,
                                'Genotype' = NA)
yearF0 = fit$data$yearF0_PCV7
vaccinationPCV7 = fit$data$yearIntroduction_PCV7
vaccinationPCV13 = fit$data$yearIntroduction_PCV13

for(k in 1:nb_countries){
  ## Compute reference fitness first:
  ## pre vacc
  mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_pre_vacc)
  yearF0_k = yearF0[k]
  vacc_PCV7_k = vaccinationPCV7[k]
  vacc_PCV13_k = vaccinationPCV13[k]
  for(gref in 1:(nb_genotypes-1)){
    if(gref == 1){
      for(i in 1){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref, yearF0_k:(vacc_PCV7_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    if(gref == 2){
      for(i in 1){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref, yearF0_k:(vacc_PCV7_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
  }
  for(i in 1:length_vect_fitness_pre_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = rep(0, nb_chains)-rep(0, nb_chains)
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,yearF0_k:(vacc_PCV7_k-1)]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(-1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  ## post vacc
  mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_post_vacc)
  for(gref in 1:(nb_genotypes-1)){
    if(gref == 1){
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] =  rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV13_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    if(gref == 2){
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] =  rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV13_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
  }
  for(i in 1:length_vect_fitness_post_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] =  rep(0, nb_chains)-rep(0, nb_chains)
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_PCV13_k:nb_years]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  
  for(g in 1:(nb_genotypes-1)){
    if(g==1){ ## PCV7
      ## pre vacc
      mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_pre_vacc)
      for(gref in 1:(nb_genotypes-1)){
        if(gref == 1){
          for(i in 1:length_vect_fitness_pre_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,1])-as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,yearF0_k:(vacc_PCV7_k-1)]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
        if(gref == 2){
          for(i in 1:length_vect_fitness_pre_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,1])-as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,yearF0_k:(vacc_PCV7_k-1)]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
      }
      for(i in 1:length_vect_fitness_pre_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,1])-rep(0, nb_chains)
        mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,yearF0_k:(vacc_PCV7_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
      }
      df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                                'Time' = rep(-1, each = nb_chains),
                                                                'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
      
      ## post vacc
      mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_post_vacc)
      for(gref in 1:(nb_genotypes-1)){
        if(gref == 1){
          for(i in 1:length_vect_fitness_post_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,1])-as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV7_k:nb_years]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
        if(gref == 2){
          for(i in 1:length_vect_fitness_post_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,1])-as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV7_k:nb_years]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
      }
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,1])-rep(0, nb_chains)
        mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_PCV7_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
      }
      df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                                'Time' = rep(1, each = nb_chains),
                                                                'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
    }
    if(g==2){ ## PCV13
      ## pre vacc
      mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_pre_vacc)
      for(gref in 1:(nb_genotypes-1)){
        if(gref == 1){
          for(i in 1:length_vect_fitness_pre_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,1])-as.vector(Chains$fitness_genotypes_pre_vacc_PCV7[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,yearF0_k:(vacc_PCV13_k-1)]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
        if(gref == 2){
          for(i in 1:length_vect_fitness_pre_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,1])-as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,yearF0_k:(vacc_PCV13_k-1)]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
      }
      for(i in 1:length_vect_fitness_pre_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_pre_vacc_PCV13[,1,1])-rep(0, nb_chains)
        mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,yearF0_k:(vacc_PCV13_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
      }
      df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                                'Time' = rep(-1, each = nb_chains),
                                                                'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
      
      ## post vacc
      mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_post_vacc)
      for(gref in 1:(nb_genotypes-1)){
        if(gref == 1){
          for(i in 1:length_vect_fitness_post_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,1])-as.vector(Chains$fitness_genotypes_post_vacc_PCV7[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV13_k:nb_years]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
        if(gref == 2){
          for(i in 1:length_vect_fitness_post_vacc){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,1])-as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,1])
            mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV13_k:nb_years]
            mean_freq = apply(mean_freq, MARGIN = 1, mean)
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
          }
        }
      }
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,1])-rep(0, nb_chains)
        mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_PCV13_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
      }
      df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                                                'Time' = rep(1, each = nb_chains),
                                                                'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
    }
  }
}

df_overall_fitness = df_overall_fitness[-1,] ## Remove the NA from the beginning

saveRDS(df_overall_fitness, '4_1_per_province_NVT_PCV7_PCV13_testswicth/disease/Data_plot_per_VT_22012024.rds')

################################################################################
## Set directory
################################################################################
setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')

################################################################################
## Fitness, per vaccine type, per vaccine era (1 switch in 2009)
################################################################################
pdf(file = 'Figure_fitness_estimates_NVT_PCV7_PCV13_disease_22012024.pdf', width = 9/2.54, height = 6/2.54)
df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
################################################################################
library(ggplot2)
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
  scale_color_manual(name = "Clades", values = c('royalblue', 'darkgreen', 'firebrick'))+                                                  
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
################################################################################

################################################################################
## Ouput estimates
################################################################################
res = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                 'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')]), digits = 2), collapse =  '-')),
                 'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')]), digits = 2), collapse =  '-'),
                                         paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')]), digits = 2), collapse =  '-')))
res_mean = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                      'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')])[1], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')])[1], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')])[1], digits = 2), collapse =  '-')),
                      'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')])[1], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')])[1], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')])[1], digits = 2), collapse =  '-')))
res_cimin = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                      'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')])[2], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')])[2], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')])[2], digits = 2), collapse =  '-')),
                      'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')])[2], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')])[2], digits = 2), collapse =  '-'),
                                              paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')])[2], digits = 2), collapse =  '-')))
res_cimax = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                       'Fitness_2000_2009' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')])[3], digits = 2), collapse =  '-'),
                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')])[3], digits = 2), collapse =  '-'),
                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')])[3], digits = 2), collapse =  '-')),
                       'Fitness_2010-2015' = c(paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')])[3], digits = 2), collapse =  '-'),
                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')])[3], digits = 2), collapse =  '-'),
                                               paste0(round(mean.and.ci(df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')])[3], digits = 2), collapse =  '-')))
saveRDS(list(res_mean, res_cimin, res_cimax), 'Fitness_estimates_NVT_PCV7_PCV13_swicth2009_plus0.rds')
################################################################################


################################################################################
## Compute VT fitness estimate - PAPER ESTIMATES
################################################################################
fitness_VT_pre = rep(0, length(nb_chains))
fitness_VT_post = rep(0, length(nb_chains))
mean_freq_PCV7_pre = NULL
mean_freq_PCV7_post = NULL
mean_freq_PCV13_pre = NULL
mean_freq_PCV13_post = NULL
for(k in 1:nb_countries){
  mean_freq = Chains$pred_absolute_freq[,k,1,vacc_WCV_k:(yearF0_k-1)]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV7_pre = c(mean_freq_PCV7_pre, mean_freq)
  
  mean_freq = Chains$pred_absolute_freq[,k,1,yearF0_k:nb_years]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV7_post = c(mean_freq_PCV7_post, mean_freq)

  mean_freq = Chains$pred_absolute_freq[,k,2,vacc_WCV_k:(yearF0_k-1)]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV13_pre = c(mean_freq_PCV13_pre, mean_freq)
  
  mean_freq = Chains$pred_absolute_freq[,k,2,yearF0_k:nb_years]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV13_post = c(mean_freq_PCV13_post, mean_freq)
}
fitness_VT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre) +
  df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre)
fitness_VT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_post/(mean_freq_PCV7_post+mean_freq_PCV13_post) +
  df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_post/(mean_freq_PCV7_post+mean_freq_PCV13_post)
fitness_NVT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')]
fitness_NVT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')]

mean.and.ci(fitness_VT_post)
mean.and.ci(fitness_NVT_post)
mean.and.ci(exp(log(fitness_VT_post)*35/365))
mean.and.ci(exp(log(fitness_NVT_post)*35/365))
mean.and.ci(fitness_NVT_post/fitness_VT_post)


mean.and.ci(fitness_VT_post/fitness_VT_pre)
mean.and.ci(fitness_NVT_post/fitness_NVT_pre)


