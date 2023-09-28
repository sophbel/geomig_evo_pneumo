## Analysis fitness Strep Pneumo serotypes by vaccine era

## Set wd 
setwd('./Fitness/')

## Load data
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
## Prepare data: each serotype, over all provinces
################################################################################
nb_clades = length(levels(as.factor(GPS_SA.sub$In_Silico_Serotype)))
nb_years = max(GPS_SA.sub$Col_Year)-min(GPS_SA.sub$Col_Year)+1
serotypes = levels(as.factor(GPS_SA.sub$In_Silico_Serotype))
nb_countries = 1
first_year = min(GPS_SA.sub$Col_Year)

threshold = 0
ref_clade = 11 ## Ref clade: Serotype 13: NVT, and has no 0s
subsample = NA
################################################################################

################################################################################
## Size arrays
################################################################################
clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_ref = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
clade_number_array_freq_non_zeros = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_non_zeros_ref = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
non_zero_country_year = matrix(0, nb_years, nb_countries)
non_zero_country_year_genotype = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
clade_number_array_freq_first_non_zeros = rep(0, nb_countries)
total_number = matrix(0, nrow = nb_countries, ncol = nb_years)
index_paramater = rep(NA, nb_clades-1)
################################################################################

################################################################################
## Fill with Strep Pneumo data
################################################################################
for(i in 1:nb_clades){
  for(j in 1:nb_years){
    clade_number_array[i,j,] = length(which(GPS_SA.sub$Col_Year == first_year+(j-1) & GPS_SA.sub$In_Silico_Serotype == serotypes[i]))
  }
}

## Compute freq with respect to a ref clade
for(kkk in 1:nb_countries){
  tmp = 1
  for(i in 1:nb_clades){
    if(i != ref_clade){
      clade_number_array_freq_ref[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[ref_clade,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref[which(is.na(clade_number_array_freq_ref) == T)] = 0

## Counts, per province, per year, for the ref clade
clade_number_array_ref = matrix(0, nrow = nb_years, ncol = nb_countries)
clade_number_array_ref[,1] = clade_number_array[ref_clade,,]

## Counts, per province, per year, per clade (without ref clade)
clade_number_array_paired = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
clade_number_array_paired[,,1] = clade_number_array[-ref_clade,,]

## Number of sequences per province
for(kkk in 1:nb_countries){
  total_number[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which province-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    non_zero_country_year_genotype[which(clade_number_array_paired[,i,kkk]+clade_number_array_ref[i,kkk] > 0),i,kkk] = 1
  }
}

## Which province-year are not 0s
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    if(sum(c(clade_number_array_paired[,i,kkk], clade_number_array_ref[i,kkk])) > 0){
      non_zero_country_year[i, kkk] = 1
    }
  }
}

## Year of vaccine introduction, per province
vaccine_introduction = rep(2009, nb_countries) - first_year + 1
index_paramater = 1:length(index_paramater)

## Gather data in 1 object for MCMC
data.MCMC = list(nb_genotypes = nb_clades,
                 nb_years = nb_years,
                 nb_countries = nb_countries,
                 
                 data_genotype_non_ref = clade_number_array_paired,
                 data_genotype_ref = clade_number_array_ref,
                 data_total_number = t(total_number),
                 non_zero_country_year = non_zero_country_year,
                 non_zero_country_year_genotype = non_zero_country_year_genotype,
                 number_zeros_country_year = length(which(non_zero_country_year == 0)),
                 number_zeros_country_year_genotype = length(which(non_zero_country_year_genotype == 0)),
                 
                 vaccine_introduction = vaccine_introduction,
                 index_paramater = index_paramater,
                 nb_vaccine_parameters = length(index_paramater))

################################################################################
## Save data 
################################################################################
saveRDS(data.MCMC, '2_processed_data/Serotypes/Data_model_individual_serotypes_11082022_ref_13NVT.rds')
################################################################################
