#######################################################
## Run fitness model including all vt-R vt-S nvt-R nvt-S
#######################################################
### Author: Sophie Belman
### Last modification: 11NOV22
#######################################################
setwd("./Fitness/")


## Load required packages
library(rstan)
library(doParallel)  
library(loo)
# library(gcc)
rstan_options(auto_write = TRUE)
############# Load model and compile
model.MCMC <- stan_model(file = './3_model/AMR_VaxStat/Model_fitness_genotypes_2p_AMR_VT_together.stan')
#
############# Data
data.MCMC = readRDS('./2_processed_data/AMR_VaxStat/Data_model_061622_ref_R.VT_fitness.rds')

############# Parameters vaccination #############################################

## Vaccine introduction (vector of length countries)
data.MCMC$yearF0 = rep(1,9)
  
## PCV introduction (vector of length countries)
data.MCMC$yearIntroduction = data.MCMC$vaccine_introduction




##################################################################################
## Run MCMC 
##################################################################################
## Name output
name_file = './4_run_model/AMR_VaxStat/output/Output_vaxstatAMR_NVT.R'


## Number of cores to use
no_cores = 1
# registerDoParallel(cores=no_cores)  
# cl = makeCluster(no_cores) 



## Seed
seed <- floor(runif(1, min = 1, max = 1E6))
# seed = 1
print(paste0('seed = ', seed))

## Function to draw initial frequencies
f0_init = function(nb_countries, nb_geno){
  res = matrix(0, ncol = nb_geno, nrow = nb_countries)
  res[,1] = rnorm(nb_countries, mean = 1-nb_geno*0.08, sd = 0.01)
  for(i in 2:nb_geno){
    res[,i] = rnorm(nb_countries, mean = 0.08, sd = 0.01)
  }
}
iters=100
# foreach(i = 1:3)  %dopar% {
  for(i in 1:3)  {

    
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC, 
                        show_messages = TRUE, 
                        # chains = 1, cores = 1,iter= iters, chain_id = i,
                        chains = 1, cores = 1,iter= 1000, chain_id = i,
                        control = list(adapt_delta = 0.97, max_treedepth = 13)) 
  fit = list(fit=fit_delay,
             data= data.MCMC)
  
  

  log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
  ## log_lik_1 contains NAs, which is a problem to compute the WAIC index and the LOO index
  ## Quick fix: remove the NAs:
  log_lik_1 = log_lik_1[,,-which(is.na(log_lik_1[1,1,]))]
  ## Then put log_lik_1 back in its initial shape array(nb_iterations/2, nb_chains,nb_datapoints)
  dims1=iters/2
  dims2=1
  dims3=length(which(!is.na(log_lik_1)))
  log_lik_1 = array( log_lik_1, dim = c(dims1, dims2, dims3)) 
  r_eff <- relative_eff(exp(log_lik_1), cores = 2)
  loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
  waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
  print(loo_1)
  print(waic_1)
  

  Chains=rstan::extract(fit$fit)
  saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
  saveRDS(data.MCMC, file = paste0(name_file, '_data_', i, '.rds'))
  
  
  m = monitor(fit$fit, print = F)
  fit$monitor = m
  saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
}
print('Reading 1')
fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
print('Reading 2')
fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
print('Reading 3')
fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))

print('Fit')
fit = NULL
fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))

fit$data = fit1$data

print('Chains')
Chains = rstan::extract(fit$fit)

print('Writing fit')
saveRDS(fit, paste0(name_file, '_fit_all.rds'))

print('Writing chains')
saveRDS(Chains, paste0(name_file, '_chains_all.rds'))

################################################################################
## Compute WAIC and LOO of the model
################################################################################

##### Lower values=better fits
####compare waics or loos with different vaccine introduction years.

library(loo)
log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
print(waic_1)
################################################################################

