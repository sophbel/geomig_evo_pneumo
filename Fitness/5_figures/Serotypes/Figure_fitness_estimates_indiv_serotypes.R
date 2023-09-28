## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)

setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/')

## Load fit
fit = readRDS(file = '4_2_per_serotype_overvallSA/Output_individual_serotypes_swicth2009_fit_all.rds')

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
## Growth rate (fitness) weighted by proportions of serotypes
################################################################################
load('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/1_raw_data/GPS_SA.GPSC.RData')

vaccination = fit$data$booster_introduction
nb_chains = length(Chains$lp__)
nb_countries = fit$data$nb_countries
nb_years = fit$data$nb_years
length_vect_fitness_post_vacc = fit$data$number_R_post_vacc
length_vect_fitness_pre_vacc = fit$data$number_R_pre_vacc
index_parameter = fit$data$index_paramater
nb_genotypes = fit$data$nb_genotypes
ref_genotype = 11
ref_genotype_header = '13'
Genotype = levels(as.factor(GPS_SA.sub$In_Silico_Serotype))[-ref_genotype]

df_overall_fitness = data.frame('Values' = NA,
                                'Time' = NA,
                                'Genotype' = NA)

vaccinationACV = fit$data$yearIntroduction
vaccinationWCV = fit$data$yearF0

for(k in 1:nb_countries){
  ## Compute reference fitness first:
  ## pre vacc
  mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_pre_vacc)
  vacc_ACV_k = vaccinationACV[k]
  vacc_WCV_k = vaccinationWCV[k]
  if(vacc_ACV_k>=nb_years) vacc_ACV_k = nb_years-1;
  if(vacc_WCV_k<=0) vacc_WCV_k = 1;
  if(vacc_WCV_k == (vacc_ACV_k-1)) vacc_WCV_k = vacc_ACV_k-2;
  if(vacc_WCV_k >= (vacc_ACV_k)) vacc_WCV_k = vacc_ACV_k-2;
  for(gref in 1:(nb_genotypes-1)){
    for(i in 1:length_vect_fitness_pre_vacc){
      if(index_parameter[gref] > 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[gref],i])
      }
      if(index_parameter[gref] == 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-rep(0, nb_chains)
      }
      # mean_freq = Chains$pred_absolute_freq[,k,gref,(vaccination[k]-((length_vect_fitness_pre_vacc-i+1)*fit$data$R_every_pre_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_pre_vacc))]
      mean_freq = Chains$pred_absolute_freq[,k,gref, vacc_WCV_k:(vacc_ACV_k-1)]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
    }
  }
  for(i in 1:length_vect_fitness_pre_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = rep(0, nb_chains)-rep(0, nb_chains)
    # mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,(vaccination[k]-((length_vect_fitness_pre_vacc-i+1)*fit$data$R_every_pre_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_pre_vacc))]
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_WCV_k:(vacc_ACV_k-1)]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(-1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  ## post vacc
  mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_post_vacc)
  for(gref in 1:(nb_genotypes-1)){
    for(i in 1:length_vect_fitness_post_vacc){
      if(index_parameter[gref] > 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[gref],i])
      }
      if(index_parameter[gref] == 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-rep(0, nb_chains)
      }
      # mean_freq = Chains$pred_absolute_freq[,k,gref,(vaccination[k]-((length_vect_fitness_post_vacc-i+1)*fit$data$R_every_post_vacc)) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_post_vacc)-1)]
      mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_ACV_k:nb_years]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
    }
  }
  for(i in 1:length_vect_fitness_post_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] =  rep(0, nb_chains)-rep(0, nb_chains)
    # mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,(vaccination[k]-((length_vect_fitness_post_vacc-i+1)*fit$data$R_every_post_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_post_vacc))]
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_ACV_k:nb_years]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  for(g in 1:(nb_genotypes-1)){
    print(g)
    ## pre vacc
    mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_pre_vacc)
    for(gref in 1:(nb_genotypes-1)){
      for(i in 1:length_vect_fitness_pre_vacc){
        if(index_parameter[g] > 0){
          if(index_parameter[gref] > 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[g],i])-as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[gref],i])
          }
          if(index_parameter[gref] == 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[g],i])-rep(0, nb_chains)
          }
        }
        if(index_parameter[g] == 0){
          if(index_parameter[gref] > 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[gref],i])
          }
          if(index_parameter[gref] == 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-rep(0, nb_chains)
          }
        }
        # mean_freq = Chains$pred_absolute_freq[,k,gref,(vaccination[k]-((length_vect_fitness_pre_vacc-i+1)*fit$data$R_every_pre_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_pre_vacc))]
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_WCV_k:(vacc_ACV_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    for(i in 1:length_vect_fitness_pre_vacc){
      if(index_parameter[g] > 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_pre_vacc[,index_parameter[g],i])-rep(0, nb_chains)
      }
      if(index_parameter[g] == 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = rep(0, nb_chains)-rep(0, nb_chains)
      }
      # mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,(vaccination[k]-((length_vect_fitness_pre_vacc-i+1)*fit$data$R_every_pre_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_pre_vacc))]
      mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_WCV_k:(vacc_ACV_k-1)]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
    }
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                              'Time' = rep(-1, each = nb_chains),
                                              'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
    
    ## post vacc
    mat_tmp = matrix(NA, ncol =  fit$data$nb_genotypes, nrow = nb_chains*length_vect_fitness_post_vacc)
    for(gref in 1:(nb_genotypes-1)){
      for(i in 1:length_vect_fitness_post_vacc){
        if(index_parameter[g] > 0){
          if(index_parameter[gref] > 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[g],i])-as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[gref],i])
          }
          if(index_parameter[gref] == 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[g],i])-rep(0, nb_chains)
          }
        }
        if(index_parameter[g] == 0){
          if(index_parameter[gref] > 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[gref],i])
          }
          if(index_parameter[gref] == 0){
            mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-rep(0, nb_chains)
          }
        }
        # mean_freq = Chains$pred_absolute_freq[,k,gref,(vaccination[k]-((length_vect_fitness_post_vacc-i+1)*fit$data$R_every_post_vacc)) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_post_vacc)-1)]
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_ACV_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    for(i in 1:length_vect_fitness_post_vacc){
      if(index_parameter[g] > 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_post_vacc[,index_parameter[g],i])-rep(0, nb_chains)
      }
      if(index_parameter[g] == 0){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = rep(0, nb_chains)-rep(0, nb_chains)
      }
      # mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,(vaccination[k]-((length_vect_fitness_post_vacc-i+1)*fit$data$R_every_post_vacc)+1) : (vaccination[k]-((length_vect_fitness_pre_vacc-i)*fit$data$R_every_post_vacc))]
      mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_ACV_k:nb_years]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
    }
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                              'Time' = rep(1, each = nb_chains),
                                              'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
  }
}

df_overall_fitness = df_overall_fitness[-1,] ## Remove the NA from the beginning

################################################################################
## Produce dataframe with mean estimates by vaccine type
################################################################################
df_overall_fitness_NVT_VT_grouped  = data.frame('Values' = NA,
                                                'Time' = NA,
                                                'Genotype' = NA,
                                                'Serotype' = NA)

pcv7.type <- c("4","6B", "9V","14","18C", "19F", "23F")
pcv13.type <- c("1","3","5","6A","7F","19A")
Genotype = c(levels(as.factor(GPS_SA.sub$In_Silico_Serotype))[-ref_genotype], levels(as.factor(GPS_SA.sub$In_Silico_Serotype))[ref_genotype])
nvt.types = Genotype[which(is.na(match(Genotype, c(pcv7.type, pcv13.type))))]

## PCV 7
tmp_pre_vacc = tmp_post_vacc = mean_freq_pre_vacc_tot = rep(0, nb_chains)
# Pre vacc
for(i in 1:length(pcv7.type)){
  index = which(Genotype == pcv7.type[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_WCV_k:(vacc_ACV_k-1)], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_pre_vacc = tmp_pre_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == pcv7.type[i] & df_overall_fitness$Time == -1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_pre_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(-1, each = nb_chains),
                                                                                        'Genotype' = rep('PCV7', nb_chains),
                                                                                        'Serotype' = rep('PCV7_mean', nb_chains)))
# Post vacc
mean_freq_pre_vacc_tot = rep(0, nb_chains)
for(i in 1:length(pcv7.type)){
  index = which(Genotype == pcv7.type[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_ACV_k:nb_years], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_post_vacc = tmp_post_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == pcv7.type[i] & df_overall_fitness$Time == 1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_post_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(1, each = nb_chains),
                                                                                        'Genotype' = rep('PCV7', nb_chains),
                                                                                        'Serotype' = rep('PCV7_mean', nb_chains)))

## PCV 13
tmp_pre_vacc = tmp_post_vacc = mean_freq_pre_vacc_tot = rep(0, nb_chains)
# Pre vacc
for(i in 1:length(pcv13.type)){
  index = which(Genotype == pcv13.type[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_WCV_k:(vacc_ACV_k-1)], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_pre_vacc = tmp_pre_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == pcv13.type[i] & df_overall_fitness$Time == -1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_pre_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(-1, each = nb_chains),
                                                                                        'Genotype' = rep('PCV13', nb_chains),
                                                                                        'Serotype' = rep('PCV13_mean', nb_chains)))
# Post vacc
mean_freq_pre_vacc_tot = rep(0, nb_chains)
for(i in 1:length(pcv13.type)){
  index = which(Genotype == pcv13.type[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_ACV_k:nb_years], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_post_vacc = tmp_post_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == pcv13.type[i] & df_overall_fitness$Time == 1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_post_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(1, each = nb_chains),
                                                                                        'Genotype' = rep('PCV13', nb_chains),
                                                                                        'Serotype' = rep('PCV13_mean', nb_chains)))

## NVT
tmp_pre_vacc = tmp_post_vacc = mean_freq_pre_vacc_tot = rep(0, nb_chains)
# Pre vacc
for(i in 1:length(nvt.types)){
  index = which(Genotype == nvt.types[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_WCV_k:(vacc_ACV_k-1)], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_pre_vacc = tmp_pre_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == nvt.types[i] & df_overall_fitness$Time == -1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_pre_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(-1, each = nb_chains),
                                                                                        'Genotype' = rep('NVT', nb_chains),
                                                                                        'Serotype' = rep('NVT_mean', nb_chains)))
# Post vacc
mean_freq_pre_vacc_tot = rep(0, nb_chains)
for(i in 1:length(nvt.types)){
  index = which(Genotype == nvt.types[i])
  mean_freq_pre_vacc = apply(Chains$pred_absolute_freq[,1,index, vacc_ACV_k:nb_years], MARGIN = 1, mean)
  mean_freq_pre_vacc_tot = mean_freq_pre_vacc_tot + mean_freq_pre_vacc
  tmp_post_vacc = tmp_post_vacc + log(df_overall_fitness$Values[which(df_overall_fitness$Genotype == nvt.types[i] & df_overall_fitness$Time == 1)])*mean_freq_pre_vacc
}
df_overall_fitness_NVT_VT_grouped = rbind(df_overall_fitness_NVT_VT_grouped, data.frame('Values' = exp(tmp_post_vacc/mean_freq_pre_vacc_tot),
                                                                                        'Time' = rep(1, each = nb_chains),
                                                                                        'Genotype' = rep('NVT', nb_chains),
                                                                                        'Serotype' = rep('NVT_mean', nb_chains)))
df_overall_fitness_NVT_VT_grouped = df_overall_fitness_NVT_VT_grouped[-1,] ## Remove the NA from the beginning

df_overall_fitness_to_plot = df_overall_fitness
df_overall_fitness_to_plot$Serotype = df_overall_fitness_to_plot$Genotype
for(i in Genotype){
  if(is.na(match(i, pcv7.type)) == F){
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)]# + rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)]# + rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Genotype[which(df_overall_fitness_to_plot$Genotype == i)] = 'PCV7'
  }
  if(is.na(match(i, pcv13.type)) == F){
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)]# + rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)] #+ rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Genotype[which(df_overall_fitness_to_plot$Genotype == i)] = 'PCV13'
  }
  if(is.na(match(i, nvt.types)) == F){
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == -1)]# + rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)] = df_overall_fitness_to_plot$Time[which(df_overall_fitness_to_plot$Genotype == i & df_overall_fitness_to_plot$Time == 1)] #+ rnorm(1, mean = 0, sd = 0.25)
    df_overall_fitness_to_plot$Genotype[which(df_overall_fitness_to_plot$Genotype == i)] = 'NVT'
  }
}

df_overall_fitness_all = rbind(df_overall_fitness_to_plot, df_overall_fitness_NVT_VT_grouped)
################################################################################

## N
Sampling_size_by_sero = c(rowSums(fit$data$data_genotype_non_ref[,,1]), sum(fit$data$data_genotype_ref))
names(Sampling_size_by_sero) = c(rownames(clade_sero_list[[1]])[-1], rownames(clade_sero_list[[1]])[1])

plot(Sampling_size_by_sero)
sero_to_remove = names(which(Sampling_size_by_sero < 50))

df_overall_fitness_NVT_VT_grouped = df_overall_fitness_NVT_VT_grouped[which(is.na(match(df_overall_fitness_NVT_VT_grouped$Serotype, sero_to_remove)) == T), ] 
df_overall_fitness_to_plot = df_overall_fitness_to_plot[which(is.na(match(df_overall_fitness_to_plot$Serotype, sero_to_remove)) == T), ] 

################################################################################
## Set directory
################################################################################
setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')

saveRDS(df_overall_fitness_all, 'Data_plot_per_serotype_12082022.rds')
################################################################################
## Fitness, per serotype, per vaccine era (1 switch in 2010)
################################################################################
pdf(width = 8/2.54, height = 6/2.54, file = "Figure_fitness_estimates_serotypes_NVT_PCV7_PCV13_12082022.pdf", onefile = T)
df_overall_fitness_NVT_VT_grouped$Genotype = factor(df_overall_fitness_NVT_VT_grouped$Genotype, levels = c("NVT", "PCV7", "PCV13"))
df_overall_fitness_to_plot$Genotype = factor(df_overall_fitness_to_plot$Genotype, levels = c("NVT", "PCV7", "PCV13"))
ACV = ggplot(data = df_overall_fitness_NVT_VT_grouped, aes(x = Time, y=Values, fill = Genotype)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f2, data = df_overall_fitness_to_plot, aes(x = Time, y = Values, fill = Serotype), 
               position = position_dodge2(width = 1), geom = "pointrange", size = 0.01, alpha = 0.4, colour = 'grey')+ 
  scale_fill_manual(labels=levels(df_overall_fitness_to_plot$Serotype), name = "Clades", values = c(rep(adjustcolor('grey', alpha.f = 0.4), 100),  'darkgreen', 'firebrick', 'royalblue'))+
  stat_summary(fun.data = f2, data = df_overall_fitness_NVT_VT_grouped, aes(x = Time, y=Values, colour = Genotype),
               geom = "pointrange", size = 0.1)+ 
  scale_colour_manual(labels=levels(df_overall_fitness_NVT_VT_grouped$Genotype), name = "Clades", values = c('darkgreen', 'firebrick', 'royalblue'))+
  # geom_pointrange(fill='blue', color='grey', shape=21, fatten = 20, size = 5)+
  theme_classic()+
  facet_grid(.~Genotype)+
  # geom_hline(data = simple_model_R, aes(yintercept = Z))+
  # geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  # geom_vline(xintercept = -5, linetype = "longdash", color = 'grey')+
  scale_x_continuous(limits = c(-2,2), 
                     breaks = c(-1,1),
                     labels = c('2000-2009', '2010-2015'))+
  scale_y_continuous(trans = 'log', limits = c(0.5,2),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  # scale_color_manual(labels=levels(df_overall_fitness_NVT_VT_grouped$Serotype), name = "Clades", values = c(rep('grey', 55),  'darkgreen', 'firebrick', 'royalblue'))+
  # scale_fill_manual(labels=unique(df_overall_fitness_NVT_VT_grouped$Genotype), name = "Clades", values = c('darkgreen', 'firebrick', 'royalblue'))+
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 10, hjust=1, angle = 30),
        axis.text.y = element_text(size = 10, hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
ACV
dev.off()

################################################################################
## Fitness, per serotype
################################################################################
## Define pcv 7 and pcv 13 serotype
pcv7.type <- c("4","6B", "9V","14","18C","23F","19F")
pcv13.type <- c("1","3","5","6A","7F","19A")
## Add label NVT, PCV7, PCV13 
Genotype = levels(as.factor(GPS_SA.sub$In_Silico_Serotype))
titles_vt = rep('NVT', length(Genotype))
titles_vt[which(is.na(match(Genotype, pcv7.type)) == F)] = 'PCV7'
titles_vt[which(is.na(match(Genotype, pcv13.type)) == F)] = 'PCV13'
colors = rep('NVT', length(Genotype))
colors[which(titles_vt == 'NVT')] = 'royalblue'
colors[which(titles_vt == 'PCV7')] = 'darkgreen'
colors[which(titles_vt == 'PCV13')] = 'firebrick'

pdf(file = 'Figure_estimates_perserotypes_12082022.pdf', width = 19/2.54, height = 27/2.54)
ACV = ggplot(df_overall_fitness, aes(x = Time, y=Values, color = Genotype)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f2, lwd = 0.2, alpha=1, position = position_dodge2(width = 0.2), geom = "pointrange")+ 
  theme_classic()+
  facet_wrap(.~Genotype, ncol = 11)+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  scale_x_continuous(limits = c(-2,2), 
                     breaks = c(-1,1),
                     labels = c('2000-2009', '2010-2015'))+
  scale_y_continuous(trans = 'log', limits = c(0.6,1.5),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "Clades", values = colors)+
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
## Fitness, per serotype, NVT combined
################################################################################
pdf(file = 'Figure_estimates_perserotypes_combinedNVT.pdf', width = 18/2.54, height = 6/2.54)

pcv7.type <- c("4","6B", "9V","14","18C", "19F", "23F")
pcv13.type <- c("1","3","5","6A","7F","19A")
df_overall_fitness_NVT_grouped = df_overall_fitness
df_overall_fitness_NVT_grouped$Genotype[which(is.na(match(df_overall_fitness_NVT_grouped$Genotype, c(pcv7.type, pcv13.type))) == T)] = 'NVT'
df_overall_fitness_NVT_grouped$Genotype = factor(df_overall_fitness_NVT_grouped$Genotype, levels = c('NVT', pcv7.type, pcv13.type))

df_overall_fitness_NVT_VT_grouped = df_overall_fitness
df_overall_fitness_NVT_VT_grouped$Genotype[which(is.na(match(df_overall_fitness_NVT_VT_grouped$Genotype, c(pcv7.type, pcv13.type))) == T)] = 'NVT'
df_overall_fitness_NVT_VT_grouped$Genotype[which(is.na(match(df_overall_fitness_NVT_VT_grouped$Genotype, c(pcv7.type))) == F)] = 'PCV7'
df_overall_fitness_NVT_VT_grouped$Genotype[which(is.na(match(df_overall_fitness_NVT_VT_grouped$Genotype, c(pcv13.type))) == F)] = 'PCV13'
# df_overall_fitness$Genotype = factor(df_overall_fitness$Genotype, levels = c('NVT', 'PCV7', 'PCV13'))
################################################################################
ACV = ggplot(df_overall_fitness_NVT_grouped, aes(x = Time, y=Values, color = Genotype)) + 
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
  scale_y_continuous(trans = 'log', limits = c(0.6,1.5),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Fitness')+
  scale_color_manual(name = "Clades", values = c('royalblue', rep('darkgreen', 7), rep('firebrick', 6)))+
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



