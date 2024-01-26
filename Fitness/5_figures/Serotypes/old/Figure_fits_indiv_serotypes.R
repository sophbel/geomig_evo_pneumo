## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)

setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL')

## Load fit
fit = readRDS(file = '4_2_per_serotype_overvallSA/Output_individual_serotypes_swicth2009_fit_all.rds')

## Chains
Chains=rstan::extract(fit$fit)

nb_countries = fit$data$nb_countries
nb_genotypes = fit$data$nb_genotypes
nb_years = fit$data$nb_years
nb_countries = fit$data$nb_countries

ref_clade = 1 ## or 5
numbers_obs_all_clades = array(NA, dim = c(nb_genotypes, nb_years, nb_countries)) 
for(i in 1:nb_countries){
  numbers_obs_all_clades[1:(nb_genotypes-1),,i] = fit$data$data_genotype_non_ref[,,i]
  numbers_obs_all_clades[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]
}
threshold = 5

min_date = 2000
nb_chains = length(Chains$lp__)

## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

############################################################################################
## Plot options
############################################################################################
setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6

############################################################################################
## Plot fits for model with 1 switch 2007 - full parametrization
############################################################################################
## PDF
## Load data for names serotypes
load('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/1_raw_data/GPS_SA.GPSC.RData')

serotypes = levels(as.factor(GPS_SA.sub$In_Silico_Serotype))

pdf(width = 19/2.54, height = 27/2.54, file = "Figure_fits_per_serotype_full_parameters.pdf", onefile = T)
par(mfrow = c(8,6), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))

cexmain = 0.9

titles = c(serotypes[-11], serotypes[11])
colors = rep(NA, length(titles))

## Define pcv 7 and pcv 13 serotype
pcv7.type <- c("4","6B", "9V","14","18C","23F","19F")
pcv13.type <- c("1","3","5","6A","7F","19A")

## Add label NVT, PCV7, PCV13 
titles_vt = rep('NVT', length(titles))
titles_vt[which(is.na(match(titles, pcv7.type)) == F)] = 'PCV7'
titles_vt[which(is.na(match(titles, pcv13.type)) == F)] = 'PCV13'

colors[which(titles_vt == 'NVT')] = 'royalblue'
colors[which(titles_vt == 'PCV7')] = 'darkgreen'
colors[which(titles_vt == 'PCV13')] = 'firebrick'
  
for(c in 1:fit$data$nb_countries){
  count = 1
  d = fit$data$data_genotype_non_ref[count,,c]
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    t = fit$data$data_total_number[,c]
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    if(i < nb_genotypes){
      d = fit$data$data_genotype_non_ref[count,,c]
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
      d = fit$data$data_genotype_ref[,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
    }
    zeros = fit$data$non_zero_country_year[,c]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = Chains$pred_absolute_freq[,c,i,1:nb_years]
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,0.2)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,16),cex = cexdots,
           main = paste0(titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i == 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                xlab="Time", ylab = 'Proportion', ylim = ylims,xlim = c(0,16),cex = cexdots,
                                main = titles[i], cex.main = cexmain,
                                col = adjustcolor('grey30', alpha.f = 0.6), 
                                yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c < 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                              xlab="", ylab = '', ylim = ylims,xlim = c(0,16),cex = cexdots,
                              main = titles[i], cex.main = cexmain,
                              col = adjustcolor('grey30', alpha.f = 0.6), 
                              yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                               xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,16), cex = cexdots,
                               main = titles[i], cex.main = cexmain,
                               col = adjustcolor('grey30', alpha.f = 0.6), 
                               yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    
    
    arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    # axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    axis(2, las = 2, at = c(0,0.1,0.2), labels = c(0,0.1,0.2), cex.axis = cexaxis, lwd = 0.9)
    axis(1, at = c(0,5, 10, 15)+1, c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
    
    lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
          cex=1.4,xlab="",ylab="", 
          col = adjustcolor(colors[i], alpha.f = 0.7))
    polygon(x = c(1:nb_years, rev(1:nb_years)),
            y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colors[i], alpha.f = 0.4), border = F)
    
    # abline(v = fit$data$yearF0[c], lty = 3)
    abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
  # mtext(Cnty[c], outer = T, cex = 1)
}
dev.off()
############################################################################################






