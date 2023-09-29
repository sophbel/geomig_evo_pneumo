############################################################################################
## Data creator - GPSC-sero model
## Noemie Lefrancq
## Last update 28/09/2023
############################################################################################

################################################################################
## Set wd 
################################################################################
# setwd('~/Documents/Project_fitness_Strep_Pneumo/geomig_evo_pneumo/Fitness')
setwd('geomig_evo_pneumo/Fitness')

################################################################################
## Load data
################################################################################
load('1_raw_data/GPS_SA.GPSC.RData')

## Define pcv 7 and pcv 13 serotype
pcv7.type <- c("4","6B", "9V","14","18C","23F","19F")
pcv13.type <- c("1","3","5","6A","7F","19A")

## Add label NVT, PCV7, PCV13 
GPS_SA.sub$Vaccine_type = NA
GPS_SA.sub$Vaccine_type = rep('NVT', length(GPS_SA.sub$Vaccine_type))
GPS_SA.sub$Vaccine_type[which(is.na(match(GPS_SA.sub$In_Silico_Serotype, pcv7.type)) == F)] = 'PCV7'
GPS_SA.sub$Vaccine_type[which(is.na(match(GPS_SA.sub$In_Silico_Serotype, pcv13.type)) == F)] = 'PCV13'
################################################################################

################################################################################
## Prepare data: find GPSC that are present at at least 1%
################################################################################
GPSCs = levels(as.factor(GPS_SA.sub$new_GPSC))
nb_per_GPSCs = NULL
for(i in 1:length(GPSCs)){
  nb_per_GPSCs = c(nb_per_GPSCs, length(which(GPS_SA.sub$new_GPSC == GPSCs[i])))
}
nb_per_GPSCs_prop = nb_per_GPSCs/sum(nb_per_GPSCs)

## Filter GPSCs
a = which(nb_per_GPSCs_prop >= 0.01)
GPSCs_sub = GPSCs[a] ## GPSCs to use 

## Selected dataset
GPS_SA.sub.selected = GPS_SA.sub[which(is.na(match(GPS_SA.sub$new_GPSC, GPSCs_sub)) == F),]
################################################################################

################################################################################
## Create array with number of strains per GPSC-sero per year
################################################################################
years = summary(GPS_SA.sub.selected$Col_Year)[1]:summary(GPS_SA.sub.selected$Col_Year)[6]
GPSCs = levels(as.factor(GPS_SA.sub.selected$new_GPSC))
seros = levels(as.factor(GPS_SA.sub.selected$In_Silico_Serotype))

mat_GPSCs = matrix(NA, ncol = length(years), nrow = length(GPSCs)*length(seros))

for(i in 1:length(GPSCs)){
  for(j in 1:length(seros)){
    print(j+(i-1)*length(seros))
    for(k in 1:length(years)){
      mat_GPSCs[j+(i-1)*length(seros),k] = length(which(GPS_SA.sub.selected$new_GPSC == GPSCs[i] & 
                                                        GPS_SA.sub.selected$In_Silico_Serotype == seros[j] &
                                                        GPS_SA.sub.selected$Col_Year == years[k]))
    }
  }
}
names = NULL
for(i in 1:length(GPSCs)){
  names = c(names, paste0(GPSCs[i], '_', seros))
}
rownames(mat_GPSCs) = names

a = which(rowSums(mat_GPSCs)>=1) ## Minimum of 5 sequences by group of GPSCxSerotype
mat_GPSCs = mat_GPSCs[a,]
names = names[a]
################################################################################


################################################################################
## Rearrange GPSC so that the reference ones are at the end of the dataset
################################################################################
## For all models: reference = GPSC-pair "52_13" (NVT serotype)
tmp = unlist(lapply(row.names(mat_GPSCs), function(x)strsplit(x, split = '_')[[1]][1]))
mat_GPSCs = rbind(mat_GPSCs[-which(tmp == '52'),], mat_GPSCs[which(tmp == '52'),])
if(length(which(tmp == '52')) > 1) names = rownames(mat_GPSCs)
if(length(which(tmp == '52')) == 1) names = c(names[-which(tmp == '52')], names[which(tmp == '52')])
################################################################################

################################################################################
## Prepare data: each genotype / serotype
################################################################################
nb_genetoypes = length(GPSCs_sub)
nb_sero = length(seros)
dim_data = length(names)
nb_years = length(years)
nb_countries = 1
first_year = 2000

threshold = 0

ref_clade = which(names == '52_13')

subsample = NA
################################################################################

################################################################################
## Size arrays
################################################################################
clade_number_array = array(0, dim = c(dim_data, nb_years))
clade_number_array_freq = array(0, dim = c(dim_data, nb_years))
clade_number_array_freq_ref = array(0, dim = c(dim_data-1, nb_years))
clade_number_array_freq_non_zeros = array(0, dim = c(dim_data, nb_years))
clade_number_array_freq_non_zeros_ref = array(0, dim = c(dim_data-1, nb_years))
non_zero_country_year = rep(0, nb_years)
non_zero_country_year_genotype = array(0, dim = c(dim_data-1, nb_years))
total_number = rep(0, nb_years)
################################################################################

################################################################################
## Fill with Strep Pneumo data
################################################################################
clade_number_array = mat_GPSCs

## Compute freq with respect to a ref clade
tmp = 1
for(i in 1:(dim_data)){
  if(i != ref_clade){
    clade_number_array_freq_ref[tmp,] = clade_number_array[i,]/(clade_number_array[ref_clade,]+clade_number_array[i,]);
    tmp = tmp+1
  }
}
clade_number_array_freq_ref[which(is.na(clade_number_array_freq_ref) == T)] = 0

## Counts, per province, per year, for the ref clade
clade_number_array_ref = clade_number_array[ref_clade,]
## Counts, per province, per year, per clade (without ref clade)
clade_number_array_paired = clade_number_array[-ref_clade,]

## Number of sequences per province
total_number =  as.vector(apply(clade_number_array, MARGIN = 2, sum))

## Which province-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(i in 1:nb_years){
  non_zero_country_year_genotype[which(clade_number_array_paired[,i]+clade_number_array_ref[i] > 0),i] = 1
}

## Which province-year are not 0s
for(i in 1:nb_years){
  if(sum(c(clade_number_array_paired[,i], clade_number_array_ref[i])) > 0){
    non_zero_country_year[i] = 1
  }
}

## Year of vaccine introduction, per province
vaccine_introduction = 2009 + 1 - first_year


## Fitness from previous fit (un-comment the relevant part)

## SEROTYPE fitness from previous fit
## Per serotype fitness, with reference serotype 13
# fit = readRDS(file = '4_run_model/Serotypes/Output_individual_serotypes_swicth2009_fit_all.rds')
# Chains=rstan::extract(fit$fit)
# sero_names_data = unlist(lapply(names[-ref_clade], function(x)stringr::str_split(x, pattern = '_')[[1]][2]))
# sero_names_fit = levels(as.factor(GPS_SA.sub$In_Silico_Serotype))
# sero_names_fit = sero_names_fit[-11]
# a = match(sero_names_data, sero_names_fit)
# fitness_genotypes_vector_pre_vacc = apply(Chains$fitness_genotypes_pre_vacc[,a,1], MARGIN = 2, mean)
# fitness_genotypes_vector_post_vacc = apply(Chains$fitness_genotypes_post_vacc[,a,1], MARGIN = 2, mean)

## VT fitness from previous fit
## Per serotype fitness, with reference serotype NVT (serotype 13 is NVT)
# fit = readRDS(file = '4_run_model/NVT_PCV7_PCV13/Output_per_provice_NVT_PCV7_PCV13_swicth2009_plus0_fit_all.rds')
# Chains=rstan::extract(fit$fit)
# sero_names = unlist(lapply(names[-ref_clade], function(x)stringr::str_split(x, pattern = '_')[[1]][2]))
# sero_names[which(is.na(match(sero_names, pcv13.type)) == T & is.na(match(sero_names, pcv7.type)) == T)] = 'NVT'
# sero_names[which(is.na(match(sero_names, pcv7.type)) == F)] = 'PCV7'
# sero_names[which(is.na(match(sero_names, pcv13.type)) == F)] = 'PCV13'
# a = match(sero_names, c('PCV7', 'PCV13')) ## PCV7 PCV13 is the order in the estimated fitness vector
# fitness_genotypes_vector_pre_vacc = apply(Chains$fitness_genotypes_pre_vacc[,a,1], MARGIN = 2, mean)
# fitness_genotypes_vector_post_vacc = apply(Chains$fitness_genotypes_vector_post_vacc[,a,1], MARGIN = 2, mean)
# fitness_genotypes_vector_pre_vacc[which(is.na(fitness_genotypes_vector_pre_vacc))] = 50 ## To mark reference
# fitness_genotypes_vector_post_vacc[which(is.na(fitness_genotypes_vector_post_vacc))] = 50 ## To mark reference

## No fitness: everything to 0
fitness_genotypes_vector_pre_vacc = rep(0, length(names)-1)
fitness_genotypes_vector_post_vacc = rep(0, length(names)-1)

## order for the covariates
data.MCMC = list(nb_GPSC = nb_genetoypes,
                 nb_sero = nb_sero,
                 dim_data = dim_data,
                 
                 nb_years = nb_years,
                 nb_countries = nb_countries,
                 
                 data_genotype_non_ref = clade_number_array_paired,
                 data_genotype_ref = clade_number_array_ref,
                 data_total_number = total_number,
                 non_zero_country_year = non_zero_country_year,
                 non_zero_country_year_genotype = non_zero_country_year_genotype,
                 number_zeros_country_year = length(which(non_zero_country_year == 0)),
                 number_zeros_country_year_genotype = length(which(non_zero_country_year_genotype == 0)),
                 
                 vaccine_introduction = vaccine_introduction, 
                 
                 fitness_genotypes_post_vacc = fitness_genotypes_vector_post_vacc,
                 fitness_genotypes_pre_vacc = fitness_genotypes_vector_pre_vacc,  
                 
                 GPSC_names = c(names[-ref_clade], names[ref_clade]))

################################################################################

################################################################################
## Write output
################################################################################
setwd('2_processed_data/GPSC-serotypes/')

## For sero fitness
# saveRDS(data.MCMC, 'Data_model_GPSC-sero_12092023_ref_NVT_GPSC_52_13_inputfitnessSERO.rds')

## For VT fitness
# saveRDS(data.MCMC, 'Data_model_GPSC-sero_12092023_ref_NVT_GPSC_52_13_inputfitnessVT.rds')

## For zero fitness
saveRDS(data.MCMC, 'Data_model_GPSC-sero_12092023_ref_NVT_GPSC_52_13_inputfitnessZERO.rds')
################################################################################

