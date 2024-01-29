## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(stringr)

setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/')

############################################################################################
## Load data
############################################################################################
## Load fit
fit_zero = readRDS(file = '4_3_per_GPSc-sero_pair_using_fitness_estimated/Fit_with_zero_fitness/Output_zero_fitness_swicth2009_fit_all.rds')
fit_VT = readRDS(file = '4_3_per_GPSc-sero_pair_using_fitness_estimated/Fit_with_VT_fitness/Output_per_VT_swicth2009_fit_all.rds')
fit_sero = readRDS(file = '4_3_per_GPSc-sero_pair_using_fitness_estimated/Fit_with_sero_fitness/Output_per_serotype_swicth2009_fit_all.rds')

## Chains
Chains_zero=rstan::extract(fit_zero$fit)
Chains_VT=rstan::extract(fit_VT$fit)
Chains_sero=rstan::extract(fit_sero$fit)
############################################################################################

############################################################################################
## Shared parameters
############################################################################################
nb_genotypes = fit_VT$data$dim_data
dim_data = fit_VT$data$dim_data
nb_years = fit_VT$data$nb_years
nb_countries = fit_VT$data$nb_countries
nb_GPSC = fit_VT$data$nb_GPSC

threshold = 5

min_date = 2000
nb_chains = length(Chains$lp__)

## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}
############################################################################################

############################################################################################
## Compute AICs
############################################################################################
AIC_table = data.frame('Model' = c('No fitness', 'VT', 'Serotype'),
                        'Number parameters' = c(((fit_zero$data$dim_data-1)), 
                                                        ((fit_VT$data$dim_data-1) + 2*(length(unique(fit_VT$data$fitness_genotypes_post_vacc))-1)), 
                                                        ((fit_sero$data$dim_data-1) + 2*(length(unique(fit_sero$data$fitness_genotypes_post_vacc))-1))),
                        'AIC' = rep(NA, 3))
AIC_table$AIC[1] = 2*((fit_zero$data$dim_data-1)) - 2 * mean(rowSums(Chains_zero$log_lik))
AIC_table$AIC[2] = 2*((fit_VT$data$dim_data-1) + 2*(length(unique(fit_VT$data$fitness_genotypes_post_vacc))-1)) - 2 *mean(rowSums(Chains_VT$log_lik))
AIC_table$AIC[3] = 2*((fit_sero$data$dim_data-1) + 2*(length(unique(fit_sero$data$fitness_genotypes_post_vacc))-1)) - 2 *mean(rowSums(Chains_sero$log_lik))
############################################################################################

############################################################################################
## Plot options
############################################################################################
setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 1

colour_sero = 'chocolate3'
colour_VT = 'firebrick'
colour_nofitness = 'chartreuse4'
############################################################################################

############################################################################################
## Plot GPSC proportions, sero model
############################################################################################
pdf(width = 25/2.54, height = 20/2.54, file = "Figure_GPSCfreq_fitness_by_serotype_23012024.pdf", onefile = T)
par(mfrow = c(5,6), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_sero
fit = fit_sero
titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  
  d_m = rep(0, nb_years)
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  for(j in 1:nb_years){
    if(t[j]>0) {
      tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
      d_m[j] = tmp$mean
      d_ci[1,j] = tmp$lower
      d_ci[2,j] = tmp$upper
    }
  }
  
  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.25)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles_GPSC_tmp_unique[tt], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_sero, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_sero, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by each serotype",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot GPSC proportions, VT model
############################################################################################
pdf(width = 25/2.54, height = 20/2.54, file = "Figure_GPSCfreq_fitness_by_VT_23012024.pdf", onefile = T)
par(mfrow = c(5,6), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_VT
fit = fit_VT
titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  
  d_m = rep(0, nb_years)
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  for(j in 1:nb_years){
    if(t[j]>0) {
      tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
      d_m[j] = tmp$mean
      d_ci[1,j] = tmp$lower
      d_ci[2,j] = tmp$upper
    }
  }
  
  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.25)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles_GPSC_tmp_unique[tt], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_VT, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_VT, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by Vaccine-Type",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot GPSC proportions, no fitness model
############################################################################################
pdf(width = 25/2.54, height = 20/2.54, file = "Figure_GPSCfreq_no_fitness_23012024.pdf", onefile = T)
par(mfrow = c(5,6), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_zero
fit = fit_zero
titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  
  d_m = rep(0, nb_years)
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  for(j in 1:nb_years){
    if(t[j]>0) {
      tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
      d_m[j] = tmp$mean
      d_ci[1,j] = tmp$lower
      d_ci[2,j] = tmp$upper
    }
  }

  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.25)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles_GPSC_tmp_unique[tt], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_nofitness, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_nofitness, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("No fitness",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot actual fits, VT model
############################################################################################
pdf(width = 40/2.54, height = 40/2.54, file = "Figure_fits_fitness_by_VT_23012024.pdf", onefile = T)
par(mfrow = c(10,11), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_VT
fit = fit_VT
titles = fit$data$GPSC_names

count = 1
d = fit$data$data_genotype_non_ref[count,]
total_m = rep(0, length(d))
total_cimin = rep(0, length(d))
total_cimax = rep(0, length(d))
for(i in 1:(nb_genotypes)){
  t = fit$data$data_total_number
  d_m = rep(0, length(d))
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  if(i < nb_genotypes){
    d = fit$data$data_genotype_non_ref[count,]
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
  }
  if(i == nb_genotypes){
    d = fit$data$data_genotype_ref
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
  }
  zeros = fit$data$non_zero_country_year
  
  d_m[which(zeros==0)] = NA
  d_ci[1,which(zeros==0)] = NA
  d_ci[2,which(zeros==0)] = NA
  
  f = Chains$pred_absolute_freq[,i,1:nb_years]
  
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.2)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles[i], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  
  # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_VT, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_VT, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by Vaccine-Type",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot actual fits, serotype model
############################################################################################
pdf(width = 40/2.54, height = 40/2.54, file = "Figure_fits_fitness_by_serotype_23012024.pdf", onefile = T)
par(mfrow = c(10,11), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_sero
fit = fit_sero
titles = fit$data$GPSC_names

count = 1
d = fit$data$data_genotype_non_ref[count,]
total_m = rep(0, length(d))
total_cimin = rep(0, length(d))
total_cimax = rep(0, length(d))
for(i in 1:(nb_genotypes)){
  t = fit$data$data_total_number
  d_m = rep(0, length(d))
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  if(i < nb_genotypes){
    d = fit$data$data_genotype_non_ref[count,]
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
  }
  if(i == nb_genotypes){
    d = fit$data$data_genotype_ref
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
  }
  zeros = fit$data$non_zero_country_year
  
  d_m[which(zeros==0)] = NA
  d_ci[1,which(zeros==0)] = NA
  d_ci[2,which(zeros==0)] = NA
  
  f = Chains$pred_absolute_freq[,i,1:nb_years]
  
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.2)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles[i], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  
  # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_sero, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_sero, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by each serotype",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot actual fits,  no fitness model
############################################################################################
pdf(width = 40/2.54, height = 40/2.54, file = "Figure_fits_no_fitness_23012024.pdf", onefile = T)
par(mfrow = c(10,11), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
Chains = Chains_zero
fit = fit_zero
titles = fit$data$GPSC_names

count = 1
d = fit$data$data_genotype_non_ref[count,]
total_m = rep(0, length(d))
total_cimin = rep(0, length(d))
total_cimax = rep(0, length(d))
for(i in 1:(nb_genotypes)){
  t = fit$data$data_total_number
  d_m = rep(0, length(d))
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  if(i < nb_genotypes){
    d = fit$data$data_genotype_non_ref[count,]
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
  }
  if(i == nb_genotypes){
    d = fit$data$data_genotype_ref
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
  }
  zeros = fit$data$non_zero_country_year
  
  d_m[which(zeros==0)] = NA
  d_ci[1,which(zeros==0)] = NA
  d_ci[2,which(zeros==0)] = NA
  
  f = Chains$pred_absolute_freq[,i,1:nb_years]
  
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.2)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles[i], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  
  # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_nofitness, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_nofitness, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("No fitness",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot Observed vs expected GPSCs
############################################################################################
pdf(width = 15/2.54, height = 6/2.54, file = "Figure_expected_vs_observed_GPSCprop_23012024.pdf", onefile = T)
par(mfrow = c(1,3), oma = c(1,1,1,0), mai = c(0.5,0.5,0.4,0.1))
fit = fit_zero
Chains = Chains_zero

titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

obs = pred = matrix(NA, nrow = length(titles_GPSC_tmp_unique), ncol = nb_years)
for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  obs[tt,] = d/t

  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  pred[tt,] = apply(f, MARGIN = 2, function(x)mean(x))
}
pred[which(obs == 0)] = NA

R2_zero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('No fitness \n R2 = ', R2_zero),
     col = colour_nofitness)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_VT
Chains = Chains_VT

titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

obs = pred = matrix(NA, nrow = length(titles_GPSC_tmp_unique), ncol = nb_years)
for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  obs[tt,] = d/t
  
  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  pred[tt,] = apply(f, MARGIN = 2, function(x)mean(x))
}
pred[which(obs == 0)] = NA

R2_VT = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('VT fitness \n R2 = ', R2_VT),
     col = colour_VT)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_sero
Chains = Chains_sero

titles = fit$data$GPSC_names
titles_GPSC_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][1]))
titles_GPSC_tmp_unique = unique(titles_GPSC_tmp)
t = fit$data$data_total_number

obs = pred = matrix(NA, nrow = length(titles_GPSC_tmp_unique), ncol = nb_years)
for(tt in 1:length(titles_GPSC_tmp_unique)){
  idx = which(titles_GPSC_tmp == titles_GPSC_tmp_unique[tt])
  d = rep(0, nb_years)
  
  for(i in 1:length(idx)){
    if(idx[i] < nb_genotypes){
      d_tmp = fit$data$data_genotype_non_ref[idx[i],]
      d = d + d_tmp
    }
    if(idx[i] == nb_genotypes){
      d_tmp = fit$data$data_genotype_ref
      d = d + d_tmp
    }
  }
  obs[tt,] = d/t
  
  f = matrix(0, nrow = nb_chains, ncol = nb_years)
  for(i in 1:length(idx)){
    f_tmp = Chains$pred_absolute_freq[,idx[i],1:nb_years]
    f = f + f_tmp
  }
  f[which(is.infinite(f)==T)] = NA
  pred[tt,] = apply(f, MARGIN = 2, function(x)mean(x))
}
pred[which(obs == 0)] = NA

R2_sero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('Sero fitness \n R2 = ', R2_sero),
     col = colour_sero)
axis(2, las =2)
abline(a=0, b=1, lty = 2)
dev.off()
############################################################################################




############################################################################################
## Plot Observed vs expected OVERALL
############################################################################################
pdf(width = 15/2.54, height = 6/2.54, file = "Figure_expected_vs_observed_prop_23012024.pdf", onefile = T)
par(mfrow = c(1,3), oma = c(1,1,1,0), mai = c(0.5,0.5,0.4,0.1))
fit = fit_zero
Chains = Chains_zero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
R2_zero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('No fitness \n R2 = ', R2_zero),
     col = colour_nofitness)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_VT
Chains = Chains_VT

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
R2_VT = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('VT fitness \n R2 = ', R2_VT),
     col = colour_VT)
axis(2, las =2)
abline(a=0, b=1, lty = 2)


fit = fit_sero
Chains = Chains_sero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
R2_sero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('Sero fitness \n R2 = ', R2_sero),
     col = colour_sero)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

dev.off()
############################################################################################

pcv7.type <- c("4","6B", "9V","14","18C","23F","19F")
pcv13.type <- c("1","3","5","6A","7F","19A")

############################################################################################
## Plot Observed vs expected NVT only
############################################################################################
pdf(width = 15/2.54, height = 6/2.54, file = "Figure_expected_vs_observed_prop_NVT_23012024.pdf", onefile = T)
par(mfrow = c(1,3), oma = c(1,1,1,0), mai = c(0.5,0.5,0.4,0.1))
fit = fit_zero
Chains = Chains_zero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(is.na(match(names_tmp, c(pcv7.type, pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]

R2_zero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('No fitness \n R2 = ', R2_zero),
     col = colour_nofitness)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_VT
Chains = Chains_VT

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(is.na(match(names_tmp, c(pcv7.type, pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_VT = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('VT fitness \n R2 = ', R2_VT),
     col = colour_VT)
axis(2, las =2)
abline(a=0, b=1, lty = 2)


fit = fit_sero
Chains = Chains_sero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(is.na(match(names_tmp, c(pcv7.type, pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_sero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('Sero fitness \n R2 = ', R2_sero),
     col = colour_sero)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

dev.off()
############################################################################################

############################################################################################
## Plot Observed vs expected PCV7 only
############################################################################################
pdf(width = 15/2.54, height = 6/2.54, file = "Figure_expected_vs_observed_prop_PCV7_23012024.pdf", onefile = T)
par(mfrow = c(1,3), oma = c(1,1,1,0), mai = c(0.5,0.5,0.4,0.1))
fit = fit_zero
Chains = Chains_zero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv7.type))))

obs = obs[idx,]
pred = pred[idx,]

R2_zero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('No fitness \n R2 = ', R2_zero),
     col = colour_nofitness)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_VT
Chains = Chains_VT

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv7.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_VT = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('VT fitness \n R2 = ', R2_VT),
     col = colour_VT)
axis(2, las =2)
abline(a=0, b=1, lty = 2)


fit = fit_sero
Chains = Chains_sero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv7.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_sero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('Sero fitness \n R2 = ', R2_sero),
     col = colour_sero)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

dev.off()
############################################################################################

############################################################################################
## Plot Observed vs expected PCV13 only
############################################################################################
pdf(width = 15/2.54, height = 6/2.54, file = "Figure_expected_vs_observed_prop_PCV13_23012024.pdf", onefile = T)
par(mfrow = c(1,3), oma = c(1,1,1,0), mai = c(0.5,0.5,0.4,0.1))
fit = fit_zero
Chains = Chains_zero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]

R2_zero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('No fitness \n R2 = ', R2_zero),
     col = colour_nofitness)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

fit = fit_VT
Chains = Chains_VT

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_VT = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('VT fitness \n R2 = ', R2_VT),
     col = colour_VT)
axis(2, las =2)
abline(a=0, b=1, lty = 2)


fit = fit_sero
Chains = Chains_sero

obs = rbind(fit$data$data_genotype_non_ref, fit$data$data_genotype_ref)
t = fit$data$data_total_number
pred = obs
for(i in 1:(dim_data)){
  for(j in 1:nb_years){
    if(obs[i,j] > 0){
      pred[i,j] = mean(Chains$pred_absolute_freq[,i,j])
    }else{
      pred[i,j] = NA
    }
    obs[i,j] = obs[i,j]/t[j]
  }
}
names = rownames(obs)
names_tmp = unlist(lapply(names, function(x)str_split(x, pattern = '_')[[1]][2]))
idx = which(!is.na(match(names_tmp, c(pcv13.type))))

obs = obs[idx,]
pred = pred[idx,]
R2_sero = 1-sum((obs-pred)^2, na.rm = T)/sum((obs-mean(obs))^2, na.rm = T)
plot(obs, pred, pch = 16, bty = 'n', xlim = c(0,0.2), ylim = c(0,0.2), yaxt = 'n',
     xlab = 'Observed', ylab = 'Predicted', main = paste0('Sero fitness \n R2 = ', R2_sero),
     col = colour_sero)
axis(2, las =2)
abline(a=0, b=1, lty = 2)

dev.off()
############################################################################################

############################################################################################
## Check likelihood
############################################################################################
log_lik_R = Chains_VT$log_lik[1,]
n = 1
for(l in 1:nb_years){
  for (j in 1:(dim_data-1)){
    if(fit_VT$data$non_zero_country_year_genotype[j,l]) {
      log_lik_R[n] = 1E-30;
      log_lik_R[n] = log(dpois(fit_VT$data$data_genotype_non_ref[j,l], Chains_VT$pred_number_non_ref[1,j,l]))
      n=n+1;
    }
  }
}
log_lik_stan = Chains_VT$log_lik[1,]
plot(log_lik_R, log_lik_stan, bty = 'n',
     xlab = 'logLikelihood R',
     ylab = 'logLikelihood STAN')
abline(a = 0, b=1)
############################################################################################



R2_zero
R2_VT
R2_sero











############################################################################################
## Massive figures for supplement
############################################################################################

############################################################################################
## Plot fits frequencies of GPSC-sero (serotype order)
############################################################################################
## PDF
# grDevices::quartz(width = 20, height = 15)
pdf(width = 63/2.54, height = 16.5/2.54, file = "Figure_fits_fitness_by_serotype_orderSERO_23012023.pdf", onefile = T)
par(mfrow = c(5,22), oma = c(1,1,1,0), mai = c(0.3,0.3,0.3,0.1), mgp = c(1.1, 0.3, 0))

Chains = Chains_sero
fit = fit_sero

titles = fit$data$GPSC_names
titles_tmp = unlist(lapply(titles, function(x)str_split(x, pattern = '_')[[1]][2]))
order_titles_by_sero = order(titles_tmp)

count = 1
d = fit$data$data_genotype_non_ref[count,]
total_m = rep(0, length(d))
total_cimin = rep(0, length(d))
total_cimax = rep(0, length(d))
for(idx in 1:(nb_genotypes)){
  i = order_titles_by_sero[idx]
  t = fit$data$data_total_number
  d_m = rep(0, length(d))
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  if(i < nb_genotypes){
    d = fit$data$data_genotype_non_ref[i,]
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
  }
  if(i == nb_genotypes){
    d = fit$data$data_genotype_ref
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
  }
  zeros = fit$data$non_zero_country_year
  
  d_m[which(zeros==0)] = NA
  d_ci[1,which(zeros==0)] = NA
  d_ci[2,which(zeros==0)] = NA
  
  f = Chains$pred_absolute_freq[,i,1:nb_years]
  
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.2)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles[i], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  
  # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_sero, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_sero, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by each serotype",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

############################################################################################
## Plot fits frequencies of GPSC-sero (GPSC order)
############################################################################################
## PDF
# grDevices::quartz(width = 20, height = 15)
pdf(width = 63/2.54, height = 16.5/2.54, file = "Figure_fits_fitness_by_serotype_orderGPSC_23012023.pdf", onefile = T)
# par(mfrow = c(5,22), oma = c(3,3,3,0), mai = c(0.2,0.2,0.2,0.1))
par(mfrow = c(5,22), oma = c(1,1,1,0), mai = c(0.3,0.3,0.3,0.1), mgp = c(1.1, 0.3, 0))

Chains = Chains_sero
fit = fit_sero

titles = fit$data$GPSC_names

count = 1
d = fit$data$data_genotype_non_ref[count,]
total_m = rep(0, length(d))
total_cimin = rep(0, length(d))
total_cimax = rep(0, length(d))
for(i in 1:(nb_genotypes)){
  t = fit$data$data_total_number
  d_m = rep(0, length(d))
  d_ci = matrix(0, ncol = length(d), nrow = 2)
  if(i < nb_genotypes){
    d = fit$data$data_genotype_non_ref[i,]
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
  }
  if(i == nb_genotypes){
    d = fit$data$data_genotype_ref
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
  }
  zeros = fit$data$non_zero_country_year
  
  d_m[which(zeros==0)] = NA
  d_ci[1,which(zeros==0)] = NA
  d_ci[2,which(zeros==0)] = NA
  
  f = Chains$pred_absolute_freq[,i,1:nb_years]
  
  f[which(is.infinite(f)==T)] = NA
  f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
  
  ylims = c(0,0.2)
  plot(1:nb_years, d_m, type="p", pch=17, bty = 'n',
       xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
       main = titles[i], cex.main = cexmain,
       col = adjustcolor('grey30', alpha.f = 0.6), 
       yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
  arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor('grey30', alpha.f = 0.8))
  
  # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
  axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
  axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
  
  lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
        cex=1.4,xlab="",ylab="", 
        col = adjustcolor(colour_sero, alpha.f = 0.7))
  polygon(x = c(1:nb_years, rev(1:nb_years)),
          y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colour_sero, alpha.f = 0.4), border = F)
  
  abline(v = fit$data$yearIntroduction, lty = 3)
}
mtext("Fitness by each serotype",
      side = 3,
      line = - 0.5,
      outer = TRUE)
dev.off()
############################################################################################

