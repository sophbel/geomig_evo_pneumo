############################################################################################
## Fit mode GPSC-sero, with sero fitness
## Noemie Lefrancq
## Last update 28/09/2023
############################################################################################
library(questionr)
library(treeio)
library(ape)
library(rstan)
library(gridExtra)
library(grid)
library(gtable)
library(doParallel)  
rstan_options(auto_write = TRUE)

############################################################################################
## Set wd
############################################################################################
setwd('geomig_evo_pneumo/Fitness/4_run_model/GPSC-serotypes/Fit_with_sero_fitness/')

############# Model
model.MCMC <- stan_model(file = '../../../3_model/GPSC-serotypes/Model_fitness_GPSC-sero_pairs_2p_vaccineintro_switch_inputfitness_fit_onlprevax.stan')

############# data for MCMC ######################################################
data.MCMC = readRDS('../../../2_processed_data/GPSC-serotypes/Data_model_GPSC-sero_12092023_ref_NVT_GPSC_52_13_inputfitnessSERO.rds')

############# parameters vaccination #############################################
data.MCMC$R_every_pre_vacc = 12
data.MCMC$number_R_pre_vacc = 1; 

data.MCMC$R_every_post_vacc = 12
data.MCMC$number_R_post_vacc = 1;

## Vaccine introduction (vector of length countries)
data.MCMC$yearF0 = rep(1, data.MCMC$nb_countries)
  
## ACV introduction (vector of length countries)
data.MCMC$yearIntroduction = data.MCMC$vaccine_introduction

##################################################################################
## Run MCMC 
##################################################################################
no_cores = 3
registerDoParallel(cores=no_cores)  
cl = makeCluster(no_cores) 

name_file = 'Output_per_serotype_swicth2009'
foreach(i = 1:3)  %dopar% {
# foreach(i = 1:3)  %do% {
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC,
                        show_messages = TRUE,
                        chains = 1, cores = 1,iter= 2000, chain_id = i,
                        control = list(adapt_delta = 0.98, max_treedepth = 13))  # iter =10000 et plus de chaines  
  fit = list(fit=fit_delay,
             data= data.MCMC)
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
library(loo)
log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
print(waic_1)
################################################################################

