## Compute difference in WAIC
library(rstan)
library(loo)

#############################################################################
## Functions
#############################################################################
compute_loo = function(fitstan){
  log_lik <- extract_log_lik(fitstan$fit, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = 2)
  loo <- loo(log_lik, r_eff = r_eff, cores = 2, )
  return(loo)
}
compute_waic = function(fitstan){
  log_lik <- extract_log_lik(fitstan$fit, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = 2)
  waic <- waic(log_lik, r_eff = r_eff, cores = 2)
  return(waic)
}
#############################################################################

#############################################################################
## Filenames
#############################################################################
## Delays tested
delays = seq(-4,4,1)
filename_name = rep(NA, length(delays))

for(i in 1:length(delays)){
  if(delays[i] >= 0) name_file = paste0('Output_per_provice_NVT_PCV7_PCV13_swicth2009_plus', abs(delays[i]), '_fit_all.rds')
  if(delays[i] < 0) name_file = paste0('Output_per_provice_NVT_PCV7_PCV13_swicth2009_minus', abs(delays[i]), '_fit_all.rds')
  
  filename_name[i] = name_file
}

filename_name = c(filename_name, paste0('Output_per_provice_NVT_PCV7_PCV13_noswitch', '_fit_all.rds'))
filename_name = c(filename_name, paste0('Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011', '_fit_all.rds'))
#############################################################################

#############################################################################
## Compute WAIC and LOO
#############################################################################
list_loo = NULL
list_waic = NULL
results_loo = data.frame('N' = 1:length(filename_name),
                         'Model' = filename_name,
                         'loo' = rep(NA, length(filename_name)),
                         'p_loo' = rep(NA, length(filename_name)))

results_waic = data.frame('N' = 1:length(filename_name),
                         'Model' = filename_name,
                         'waic' = rep(NA, length(filename_name)),
                         'p_waic' = rep(NA, length(filename_name)))

for(i in 1:length(filename_name)){
  print(filename_name[i])
  
  fitstan = readRDS(file = filename_name[i])
  list_loo[[i]] = compute_loo(fitstan)
  list_waic[[i]] = compute_waic(fitstan)
  
  results_loo$Model[i] = filename_name[i]
  results_loo$loo[i] = list_loo[[i]]$estimates[3]
  results_loo$p_loo[i] = list_loo[[i]]$estimates[2]
  
  results_waic$Model[i] = filename_name[i]
  results_waic$waic[i] = list_waic[[i]]$estimates[3]
  results_waic$p_waic[i] = list_waic[[i]]$estimates[2]
}
#############################################################################

#############################################################################
## Save results
#############################################################################
saveRDS(object = list_loo, file = 'list_loo.rds')
saveRDS(object = list_waic, file = 'list_waic.rds')
write.csv(x = results_loo, file = 'results_loo.csv')
write.csv(x = results_waic, file = 'results_waic.csv')
#############################################################################