

################################################################################
###Resistance Data###
##Reference is Resistance which is the second row.
################################################################################
setwd('/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/AMR/')
load('clade_list.RData')
load('clade_numer_array.RData')
load('clade_sero_list.RData')
clade_list_amr=clade_list
# clade_number_array_strep_pneumo = clade_number_array
clade_number_array_amr = clade_number_array
clade_sero_list_amr=clade_sero_list
################################################################################
###Vaccine Type Data###
###Reference is NVT which is the 1st row
################################################################################
setwd('/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/VaccineStatus/')
load('clade_list.RData')
load('clade_numer_array.RData')
load('clade_sero_list.RData')
clade_list_vstat=clade_list
# clade_number_array_strep_pneumo = clade_number_array
clade_number_array_vstat = clade_number_array
clade_sero_list_vstat=clade_sero_list

################################################################################
###Vaccine Type Data###
###Reference is NVT-R which is the XX row
################################################################################
setwd('/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/Overall/')
load('clade_list.RData')
load('clade_numer_array.RData')
load('clade_sero_list.RData')
clade_list_total=clade_list
# clade_number_array_strep_pneumo = clade_number_array
clade_number_array_total = clade_number_array
clade_sero_list_total=clade_sero_list


################################################################################
## Prepare data: each province
################################################################################
nb_clades_vstat = dim(clade_number_array_vstat)[1]
nb_clades_amr = dim(clade_number_array_amr)[1]
nb_clades_total = dim(clade_number_array_total)[1]
nb_years = dim(clade_number_array_vstat)[2]
nb_countries = dim(clade_number_array_vstat)[3]
first_year = 2009

################################################################################
## Size arrays
################################################################################
##AMR
# clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_amr = array(0, dim = c(nb_clades_amr, nb_years, nb_countries))
clade_number_array_freq_ref_amr = array(0, dim = c(nb_clades_amr-1, nb_years, nb_countries))
clade_number_array_freq_non_zeros_amr = array(0, dim = c(nb_clades_amr, nb_years, nb_countries))
clade_number_array_freq_non_zeros_ref_amr = array(0, dim = c(nb_clades_amr-1, nb_years, nb_countries))
non_zero_country_year_amr = matrix(0, nb_years, nb_countries)
non_zero_country_year_genotype_amr = array(0, dim = c(nb_clades_amr-1, nb_years, nb_countries))
clade_number_array_freq_first_non_zeros_amr = rep(0, nb_countries)
total_number_amr = matrix(0, nrow = nb_countries, ncol = nb_years)
##VaccineStatus
# clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_vstat = array(0, dim = c(nb_clades_vstat, nb_years, nb_countries))
clade_number_array_freq_ref_vstat = array(0, dim = c(nb_clades_vstat-1, nb_years, nb_countries))
clade_number_array_freq_non_zeros_vstat = array(0, dim = c(nb_clades_vstat, nb_years, nb_countries))
clade_number_array_freq_non_zeros_ref_vstat = array(0, dim = c(nb_clades_vstat-1, nb_years, nb_countries))
non_zero_country_year_vstat = matrix(0, nb_years, nb_countries)
non_zero_country_year_genotype_vstat = array(0, dim = c(nb_clades_vstat-1, nb_years, nb_countries))
clade_number_array_freq_first_non_zeros_vstat = rep(0, nb_countries)
total_number_vstat = matrix(0, nrow = nb_countries, ncol = nb_years)
##AMR per Vaccine Status
# clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_total = array(0, dim = c(nb_clades_total, nb_years, nb_countries))
# clade_number_array_freq_ref_total = array(0, dim = c(nb_clades_total-1, nb_years, nb_countries))
clade_number_array_freq_non_zeros_total = array(0, dim = c(nb_clades_total, nb_years, nb_countries))
# clade_number_array_freq_non_zeros_ref_total = array(0, dim = c(nb_clades_total-1, nb_years, nb_countries))
non_zero_country_year_total = matrix(0, nb_years, nb_countries)
non_zero_country_year_genotype_total = array(0, dim = c(nb_clades_total, nb_years, nb_countries))
clade_number_array_freq_first_non_zeros_total = rep(0, nb_countries)
total_number_total = matrix(0, nrow = nb_countries, ncol = nb_years)
################################################################################

################################################################################
## AMR fill with Data
################################################################################
threshold = 0
ref_clade_amr = 2 ## Ref clade: Resistant
subsample = NA
clade_number_array = clade_number_array_amr

## Compute freq with respect to a ref clade
for(kkk in 1:nb_countries){
  tmp = 1
  for(i in 1:nb_clades_amr){
    if(i != ref_clade_amr){
      clade_number_array_freq_ref_amr[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[ref_clade_amr,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref_amr[which(is.na(clade_number_array_freq_ref_amr) == T)] = 0

## Counts, per province, per year, for the ref clade
clade_number_array_ref_amr = clade_number_array[ref_clade_amr,,]
## Counts, per province, per year, per clade (without ref clade)
clade_number_array_paired_amr<-array(NA,dim=c(1,nb_years,nb_countries))
for(i in 1:nb_countries){
  clade_number_array_paired_amr[,1:nb_years,i]<-clade_number_array[,,i][1,]
  
}
## Number of sequences per province
for(kkk in 1:nb_countries){
  total_number_amr[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which province-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    non_zero_country_year_genotype_amr[which(clade_number_array_paired_amr[,i,kkk]+clade_number_array_ref_amr[i,kkk] > 0),i,kkk] = 1
  }
}

## Which province-year are not 0s
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    if(sum(c(clade_number_array_paired_amr[,i,kkk], clade_number_array_ref_amr[i,kkk])) > 0){
      non_zero_country_year_amr[i, kkk] = 1
    }
  }
}
################################################################################


################################################################################
## Vaccine Status fill with Data
################################################################################
threshold = 0
ref_clade_vstat = 1 ## Ref clade: NVT
subsample = NA
clade_number_array = clade_number_array_vstat

## Compute freq with respect to a ref clade
for(kkk in 1:nb_countries){
  tmp = 1
  for(i in 1:nb_clades_vstat){
    if(i != ref_clade_vstat){
      clade_number_array_freq_ref_vstat[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[ref_clade_vstat,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref_vstat[which(is.na(clade_number_array_freq_ref_vstat) == T)] = 0

## Counts, per province, per year, for the ref clade
clade_number_array_ref_vstat = clade_number_array[ref_clade_vstat,,]
## Counts, per province, per year, per clade (without ref clade)
clade_number_array_paired_vstat<-array(NA,dim=c(1,nb_years,nb_countries))
for(i in 1:nb_countries){
  clade_number_array_paired_vstat[,1:nb_years,i]<-clade_number_array[,,i][2,]
  
}
## Number of sequences per province
for(kkk in 1:nb_countries){
  total_number_vstat[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which province-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    non_zero_country_year_genotype_vstat[which(clade_number_array_paired_vstat[,i,kkk]+clade_number_array_ref_vstat[i,kkk] > 0),i,kkk] = 1
  }
}

## Which province-year are not 0s
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    if(sum(c(clade_number_array_paired_vstat[,i,kkk], clade_number_array_ref_vstat[i,kkk])) > 0){
      non_zero_country_year_vstat[i, kkk] = 1
    }
  }
}
################################################################################

################################################################################
## AMR within Vaccine Status fill with Data
################################################################################
threshold = 0
ref_clade_total = 3 ## Ref clade: VT Resistant 
subsample = NA
clade_number_array = clade_number_array_total

## Counts, per province, per year, for the ref clade
clade_number_array_total = clade_number_array

## Number of sequences per province
for(kkk in 1:nb_countries){
  total_number_total[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which province-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    non_zero_country_year_genotype_total[which(clade_number_array_total[,i,kkk] > 0),i,kkk] = 1
  }
}

## Which province-year are not 0s
for(kkk in 1:nb_countries){
  for(i in 1:nb_years){
    if(sum(clade_number_array_total[,i,kkk]) > 0){
      non_zero_country_year_total[i, kkk] = 1
    }
  }
}
################################################################################
##proportion of each group in order S,R,NVT,VT
# freq_per_class<-matrix(nrow=nb_genotypes_total, ncol=nb_years)
#   freq_per_class[1,]<-apply(data.MCMC$data_genotype_non_ref_amr,2,sum)/rowSums(data.MCMC$data_total_number_total)
#   freq_per_class[2,]<-rowSums(data.MCMC$data_genotype_ref_amr)/rowSums(data.MCMC$data_total_number_total)
#   freq_per_class[3,]<-apply(data.MCMC$data_genotype_non_ref_vstat,2,sum)/rowSums(data.MCMC$data_total_number_total)
#   freq_per_class[4,]<-rowSums(data.MCMC$data_genotype_ref_vstat)/rowSums(data.MCMC$data_total_number_total)

  # freq_per_class<-vector(mode='numeric', length=nb_genotypes_total)
  # freq_per_class<-as.integer(c(86,100,20,166))
  # freq_per_class_all=as.integer(186)
## Year of vaccine introduction, per province
vaccine_introduction = rep(2010, nb_countries) - first_year

## order for the covariates
data.MCMC = list(nb_genotypes_amr = nb_clades_amr,
                 nb_genotypes_vstat = nb_clades_vstat,
                 nb_genotypes_total = nb_clades_total,
                 nb_years = nb_years,
                 nb_countries = nb_countries,
                 
                 data_genotype_non_ref_amr = clade_number_array_paired_amr,
                 data_genotype_ref_amr = clade_number_array_ref_amr,
                 data_total_number_amr = t(total_number_amr),
                 non_zero_country_year_amr = non_zero_country_year_amr,
                 non_zero_country_year_genotype_amr = non_zero_country_year_genotype_amr,
                 number_zeros_country_year_amr = length(which(non_zero_country_year_amr == 0)),
                 number_zeros_country_year_genotype_amr = length(which(non_zero_country_year_genotype_amr == 0)),
                 
                 data_genotype_non_ref_vstat = clade_number_array_paired_vstat,
                 data_genotype_ref_vstat = clade_number_array_ref_vstat,
                 data_total_number_vstat = t(total_number_vstat),
                 non_zero_country_year_vstat = non_zero_country_year_vstat,
                 non_zero_country_year_genotype_vstat = non_zero_country_year_genotype_vstat,
                 number_zeros_country_year_vstat = length(which(non_zero_country_year_vstat == 0)),
                 number_zeros_country_year_genotype_vstat = length(which(non_zero_country_year_genotype_vstat == 0)),
                 
                 # data_genotype_non_ref_total = clade_number_array_paired_total,
                 data_genotype_total = clade_number_array_total,
                 data_total_number_total = t(total_number_total),
                 non_zero_country_year_total = non_zero_country_year_total,
                 non_zero_country_year_genotype_total = non_zero_country_year_genotype_total,
                 number_zeros_country_year_total = length(which(non_zero_country_year_total == 0)),
                 number_zeros_country_year_genotype_total = length(which(non_zero_country_year_genotype_total == 0)),
                 # number_zeros_country_year_genotype_total = length(which(non_zero_country_year_genotype_total[-4,,] == 0)),
                 
                 # freq_per_class_all = freq_per_class_all,
                 # freq_per_class = freq_per_class,
                 
                 vaccine_introduction = vaccine_introduction)
setwd('/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/')
saveRDS(data.MCMC, 'Data_model_061622_ref_R.VT_fitness.rds')

