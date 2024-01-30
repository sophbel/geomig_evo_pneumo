## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(ggplot2)

# setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/')
## Load fit
# fit = readRDS(file = '4_1_per_province_NVT_PCV7_PCV13_testswicth/test_shift_delay_22012024/Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011_1pervax_plus0_fit_all.rds')
fit = readRDS(file='./4_run_model/NVT_PCV7_PCV13/output/Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011_1pervax_plus0_fit_all.rds')
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
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV7_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    if(gref == 2){
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] =  rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc_PCV13[,1,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_PCV7_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
  }
  for(i in 1:length_vect_fitness_post_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] =  rep(0, nb_chains)-rep(0, nb_chains)
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_PCV7_k:nb_years]
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

# saveRDS(df_overall_fitness, '4_run_model/NVT_PCV7_PCV13/output/Data_plot_per_VT_22012024.rds')
# =======
# saveRDS(df_overall_fitness, 'Data_plot_per_VT_22012024.rds')
################################################################################

################################################################################
## Compute relative fitness of each VT, compared to NVT, pre and post PCV
################################################################################
df_relative_fitness_vs_NVT = data.frame('Values' = NA,
                                        'Type' = NA,
                                        'Time' = NA)

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')
c = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('NVT', length(a)),
                                                                          'Time' = rep('Pre-PCV', length(a))))
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[b]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('PCV7', nb_chains),
                                                                          'Time' = rep('Pre-PCV', length(a))))
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[c]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('PCV13', nb_chains),
                                                                          'Time' = rep('Pre-PCV', length(a))))

a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')
c = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('NVT', length(a)),
                                                                          'Time' = rep('Post-PCV', length(a))))
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[b]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('PCV7', length(a)),
                                                                          'Time' = rep('Post-PCV', length(a))))
df_relative_fitness_vs_NVT = rbind(df_relative_fitness_vs_NVT, data.frame('Values' = df_overall_fitness$Values[c]/df_overall_fitness$Values[a],
                                                                          'Type' = rep('PCV13', length(a)),
                                                                          'Time' = rep('Post-PCV', length(a))))
df_relative_fitness_vs_NVT = df_relative_fitness_vs_NVT[-1,]
################################################################################

################################################################################
## Compute efect of vaccine implementation: relative fitness of each Type after vax implementation, compared to pre implementation
################################################################################
df_relative_fitness_type_post_vs_pre = data.frame('Values' = NA,
                                                  'Type' = NA,
                                                  'Time' = NA)

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('NVT', length(a)),
                                                                                              'Time' = rep('Pre-PCV', length(a))))
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[b]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('NVT', nb_chains),
                                                                                              'Time' = rep('Post-PCV', length(a))))

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('PCV7', length(a)),
                                                                                              'Time' = rep('Pre-PCV', length(a))))
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[b]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('PCV7', nb_chains),
                                                                                              'Time' = rep('Post-PCV', length(a))))

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('PCV13', length(a)),
                                                                                              'Time' = rep('Pre-PCV', length(a))))
df_relative_fitness_type_post_vs_pre = rbind(df_relative_fitness_type_post_vs_pre, data.frame('Values' = df_overall_fitness$Values[b]/df_overall_fitness$Values[a],
                                                                                              'Type' = rep('PCV13', nb_chains),
                                                                                              'Time' = rep('Post-PCV', length(a))))

df_relative_fitness_type_post_vs_pre = df_relative_fitness_type_post_vs_pre[-1,]
################################################################################



################################################################################
## Set directory to save figures
################################################################################
# setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')
################################################################################

################################################################################
## Fig 4D-E Relative fitness, compared to NVT, per vaccine type, per vaccine era
################################################################################
pdf(file = './5_figures/Figure4DE_relative_fitness_estimates_prevax_NVT_PCV7_PCV13_22012026.pdf', width = 4, height = 3)
df_relative_fitness_vs_NVT$Type = factor(df_relative_fitness_vs_NVT$Type, levels = c('NVT', 'PCV7', 'PCV13'))
df_relative_fitness_vs_NVT$Time = factor(df_relative_fitness_vs_NVT$Time, levels = c('Pre-PCV', 'Post-PCV'))
################################################################################
p23<-ggplot(df_relative_fitness_vs_NVT,aes(x=Type,y=Values,color=Type))+
  facet_grid(.~Time)+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  stat_summary(fun.data = f2, size=3, shape=18, position = position_dodge2(width = 0.2), geom = "point")+
  stat_summary(fun.data = f2, width=0.3, lwd = 0.5, position = position_dodge2(width = 0.2), geom = "errorbar")+ 
  theme_classic()+
  scale_y_continuous(trans = 'log', limits = c(0.5,1.6),
                     breaks = c(0.1, 0.25,0.5, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative Fitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=45,vjust=0.6),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")
p23
dev.off()
################################################################################

################################################################################
## Fig 4F Relative fitness of each Type after vax implementation, compared to pre implementation
################################################################################
pdf(file = './5_figures/Figure4F_relative_fitness_estimates_post_vs_prevax_NVT_PCV7_PCV13_22012026.pdf', width = 4, height = 3)
df_relative_fitness_type_post_vs_pre$Type = factor(df_relative_fitness_type_post_vs_pre$Type, levels = c('NVT', 'PCV7', 'PCV13'))
df_relative_fitness_type_post_vs_pre$Time = factor(df_relative_fitness_type_post_vs_pre$Time, levels = c('Pre-PCV', 'Post-PCV'))
################################################################################
p4<-ggplot(df_relative_fitness_type_post_vs_pre,aes(x=Time,y=Values,color=Type))+
  facet_grid(.~Type)+
  geom_hline(yintercept = 1, linetype = "longdash", color = 'grey')+
  stat_summary(fun.data = f2, size=3, shape=18, position = position_dodge2(width = 0.2), geom = "point")+
  stat_summary(fun.data = f2, width=0.3, lwd = 0.5, position = position_dodge2(width = 0.2), geom = "errorbar")+ 
  theme_classic()+
  scale_y_continuous(trans = 'log', limits = c(0.5,1.6),
                     breaks = c(0.1, 0.25,0.5,  1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5,  1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'RelativeFitness')+
  scale_color_manual(name = "", values = c( 'maroon','darkblue', 'darkgreen'))+
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=45,vjust=0.6),axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", angle = 0),legend.position = "none")
p4
dev.off()
################################################################################

################################################################################
## Fitness, per vaccine type, per vaccine era (NOT IN FUGURE 4!)
################################################################################
pdf(file = './5_figures/Figure_fitness_estimates_NVT_PCV7_PCV13_22012024.pdf', width = 9/2.54, height = 6/2.54)
df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
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
## Estimate relative fitness, per vax era
################################################################################
res = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                 'Pre-PCV' = c(paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Pre-PCV' & df_relative_fitness_vs_NVT$Type == 'NVT')]), digits = 2), collapse =  '-'),
                               paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Pre-PCV' & df_relative_fitness_vs_NVT$Type == 'PCV7')]), digits = 2), collapse =  '-'),
                               paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Pre-PCV' & df_relative_fitness_vs_NVT$Type == 'PCV13')]), digits = 2), collapse =  '-')),
                 'Post-PCV' = c(paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Post-PCV' & df_relative_fitness_vs_NVT$Type == 'NVT')]), digits = 2), collapse =  '-'),
                                paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Post-PCV' & df_relative_fitness_vs_NVT$Type == 'PCV7')]), digits = 2), collapse =  '-'),
                                paste0(round(mean.and.ci(df_relative_fitness_vs_NVT$Values[which(df_relative_fitness_vs_NVT$Time == 'Post-PCV' & df_relative_fitness_vs_NVT$Type == 'PCV13')]), digits = 2), collapse =  '-')))
write.csv(res, './5_figures/Estimates_Fig4DE_relative_fitness_estimates_prevax_NVT_PCV7_PCV13_22012026.csv')
################################################################################

################################################################################
## Estimate relative fitness of each Type after vax implementation, compared to pre implementation
################################################################################
res = data.frame('seroype' = c('NVT', 'PCV7', 'PCV13'),
                 'Pre-PCV' = c(paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Pre-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'NVT')]), digits = 2), collapse =  '-'),
                               paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Pre-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV7')]), digits = 2), collapse =  '-'),
                               paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Pre-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV13')]), digits = 2), collapse =  '-')),
                 'Post-PCV' = c(paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'NVT')]), digits = 2), collapse =  '-'),
                                paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV7')]), digits = 2), collapse =  '-'),
                                paste0(round(mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV13')]), digits = 2), collapse =  '-')))
write.csv(res, './5_figures/Estimates_Fig4F_relative_fitness_estimates_post_vs_prevax_NVT_PCV7_PCV13_22012026.csv')
################################################################################

################################################################################
## PAPER ESTIMATES - Compute fitness estimates for the paper
################################################################################
yearF0 = fit$data$yearF0_PCV7[k]
vaccinationPCV7 = fit$data$yearIntroduction_PCV7[k]
vaccinationPCV13 = fit$data$yearIntroduction_PCV13[k]

fitness_VT_pre = rep(0, length(nb_chains))
fitness_VT_post = rep(0, length(nb_chains))
mean_freq_PCV7_pre = NULL
mean_freq_PCV7_post = NULL
mean_freq_PCV13_pre = NULL
mean_freq_PCV13_post = NULL
for(k in 1:nb_countries){
  mean_freq = Chains$pred_absolute_freq[,k,1,yearF0:(vaccinationPCV7-1)]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV7_pre = c(mean_freq_PCV7_pre, mean_freq)
  
  mean_freq = Chains$pred_absolute_freq[,k,1,vaccinationPCV7:nb_years]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV7_post = c(mean_freq_PCV7_post, mean_freq)
  
  mean_freq = Chains$pred_absolute_freq[,k,2,yearF0:(vaccinationPCV13-1)]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV13_pre = c(mean_freq_PCV13_pre, mean_freq)
  
  mean_freq = Chains$pred_absolute_freq[,k,2,vaccinationPCV13:nb_years]
  mean_freq = apply(mean_freq, MARGIN = 1, mean)
  mean_freq_PCV13_post = c(mean_freq_PCV13_post, mean_freq)
}
fitness_VT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre) +
  df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_pre/(mean_freq_PCV7_pre+mean_freq_PCV13_pre)
fitness_VT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV7')]*mean_freq_PCV7_post/(mean_freq_PCV7_post+mean_freq_PCV13_post) +
  df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'PCV13')]*mean_freq_PCV13_post/(mean_freq_PCV7_post+mean_freq_PCV13_post)
fitness_NVT_post = df_overall_fitness$Values[which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'NVT')]
fitness_NVT_pre = df_overall_fitness$Values[which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'NVT')]

## Abstract
mean.and.ci(fitness_NVT_post/fitness_VT_post) 

## Main
mean.and.ci(fitness_NVT_pre)

a = which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV7')
mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[a]) ## OK
a = which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'PCV13')
mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[a]) ## OK
a = which(df_relative_fitness_type_post_vs_pre$Time == 'Post-PCV' & df_relative_fitness_type_post_vs_pre$Type == 'NVT')
mean.and.ci(df_relative_fitness_type_post_vs_pre$Values[a]) ## OK

mean.and.ci(fitness_NVT_post/fitness_VT_post)
mean.and.ci(exp(log(fitness_NVT_post/fitness_VT_post)*35/365))
################################################################################

