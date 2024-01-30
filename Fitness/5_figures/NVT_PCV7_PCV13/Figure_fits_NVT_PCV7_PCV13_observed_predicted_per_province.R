## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)

## Load fit
fit = readRDS(file = '4_run_model/NVT_PCV7_PCV13/output/Output_per_provice_NVT_PCV7_PCV13_swicthes2009and2011_1pervax_plus0_fit_all.rds')
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
cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6

province = c("Eastern Cape", "Free State", "Gauteng", "KwaZulu-Natal", "Limpopo", "Mpumalanga", "North West", "Northern Cape", "Western Cape")

############################################################################################
## Plot fits for model shared_param_PCV7_PCV13, per province
############################################################################################
pdf(width = 9/2.54, height = 25/2.54, file = "5_figures/NVT_PCV7_PCV13/Figure_fits_NVT_PCV7_PCV13_per_province_22012024.pdf", onefile = T)
par(mfrow = c(9,3), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
colors = c('chartreuse4', 'firebrick', 'royalblue')

for(c in 1:fit$data$nb_countries){
  titles = c('PCV7',
             'PCV13',
             'NVT')
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
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,15),cex = cexdots,
           main = paste0(province[c], '-',titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i == 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                xlab="Time", ylab = 'Proportion', ylim = ylims,xlim = c(0,15),cex = cexdots,
                                main = titles[i], cex.main = cexmain,
                                col = adjustcolor('grey30', alpha.f = 0.6), 
                                yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c < 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                              xlab="", ylab = '', ylim = ylims,xlim = c(0,15),cex = cexdots,
                              main = titles[i], cex.main = cexmain,
                              col = adjustcolor('grey30', alpha.f = 0.6), 
                              yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                               xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
                               main = titles[i], cex.main = cexmain,
                               col = adjustcolor('grey30', alpha.f = 0.6), 
                               yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    
    arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    # axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
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
}
dev.off()
############################################################################################

############################################################################################
## Plot fits for model shared_param_PCV7_PCV13 - overall south africa
############################################################################################
pdf(width = 13/2.54, height = 5/2.54, file = "5_figures/NVT_PCV7_PCV13/Figure_fits_NVT_PCV7_PCV13_overall_22012024.pdf", onefile = T)
par(mfcol = c(1,3), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))

colors = c('chartreuse4', 'firebrick', 'royalblue')
for(c in 1:1){
  titles = c('PCV7',
             'PCV13',
             'NVT')
  count = 1
  d = rowSums(fit$data$data_genotype_non_ref[count,,])
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  fvt.mat<-d_m.mat<-matrix(nrow=nb_genotypes,ncol=nb_years)
  mat.list<-list()
  
  for(i in 1:(nb_genotypes)){
    mat<-matrix(nrow=nb_years,ncol=7)
    t = rowSums(fit$data$data_total_number)
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    if(i < nb_genotypes){
      d = rowSums(fit$data$data_genotype_non_ref[count,,])
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
      d = rowSums(fit$data$data_genotype_ref)
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
    }
    zeros = fit$data$non_zero_country_year[,1]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    ##### save the proportion of each in each year
    d_m.mat[i,]<-d_m
    
    f = matrix(0, nrow = nb_chains, ncol = length(1:nb_years))
    if(i < nb_genotypes){
      for(ccc in 1:nb_countries){
        # f = f + Chains$pred_number_non_ref[,ccc,i,1:nb_years]
        # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*fit$data$data_total_number[,ccc]/t))
        # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*1/9)) 
        f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*sum(fit$data$data_total_number[,ccc])/sum(t)))
      }
    }
    if(i == nb_genotypes){
      for(ccc in 1:nb_countries){
        # f = f + Chains$pred_number_ref[,ccc,1:nb_years]
        # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*fit$data$data_total_number[,ccc]/t))
        # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*1/9))
        f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*sum(fit$data$data_total_number[,ccc])/sum(t)))
      }
      # f = t(apply(f, MARGIN = 1, function(x)x/t))
    }
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    ###store the fits across years for each type
    fvt.mat[i,]<-f_mean_ci[1,]
    
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,15),cex = cexdots,
           main = paste0(titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i == 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                xlab="Time", ylab = '   ', ylim = ylims,xlim = c(0,15),cex = cexdots,
                                main = titles[i], cex.main = cexmain,
                                col = adjustcolor('grey30', alpha.f = 0.6), 
                                yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c < 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                              xlab="", ylab = '', ylim = ylims,xlim = c(0,15),cex = cexdots,
                              main = titles[i], cex.main = cexmain,
                              col = adjustcolor('grey30', alpha.f = 0.6), 
                              yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    if(i > 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                               xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,15), cex = cexdots,
                               main = titles[i], cex.main = cexmain,
                               col = adjustcolor('grey30', alpha.f = 0.6), 
                               yaxt = 'n', xaxt = 'n', cex.lab = cexlab)}
    
    
    arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    # axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
    axis(1, at = c(0,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
    
    lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
          cex=1.4,xlab="",ylab="", 
          col = adjustcolor(colors[i], alpha.f = 0.7))
    polygon(x = c(1:nb_years, rev(1:nb_years)),
            y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colors[i], alpha.f = 0.4), border = F)
    
    # abline(v = fit$data$yearF0[c], lty = 3)
    # abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
    
    mat[,1]<-1999+1:nb_years
    mat[,2]<-t(d_m)
    mat[,3:4]<-t(d_ci)
    mat[,5:7]<-t(f_mean_ci)
    colnames(mat)<-c("years","data","datalower","dataupper","fit","fitlower","fitupper")
    mat<-data.table(mat)
    mat$type<-titles[i]
    mat.list[[i]]<-mat
  }
}
dev.off()
mat_nvtvt<-rbindlist(mat.list)
save(mat_nvtvt,file="5_figures/NVT_PCV7_PCV13/fit_table_nvtpcv7pcv13_2switch.RData")
############################################################################################



############################################################################################
## Part 2: Observed predicted
############################################################################################
rd = c(1:9) ## order plot country
cols = RColorBrewer::brewer.pal('Set1', n = 9)

cexdots= 0.8
pdf(width = 3/2.54, height = 9/2.54, file = "Figure_observed_predicted_loose_priors_22012024.pdf", onefile = T)
# windows(width = 3/2.54, height = 9/2.54)
par(mfrow = c(3,1), oma = c(0,0,0,0), mai = c(0.2,0.2,0.2,0.1))
## Numbers
numbers_obs = fit$data$data_genotype_non_ref
numbers_pred_m = numbers_pred_cimax = numbers_pred_cimin = array(0, dim = c(nb_genotypes-1, nb_years, nb_countries))
for(i in 1:nb_countries){
  for(j in 1:(nb_genotypes-1)){
    numbers_pred_m[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[1])
    numbers_pred_cimin[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[2])
    numbers_pred_cimax[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[3])
  }
  numbers_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
plot(numbers_obs[,,rd[1]], numbers_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), pch = 20, xlab = 'Observed', ylab = 'Predicted', yaxt = 'n', bty = 'n', xaxt = 'n',
     main = "Numbers",
     cex.axis = cexaxis, cex.lab = cexlab, cex = cexdots,
     xlim = c(0.5,max(Chains$pred_number_non_ref)), ylim = c(0.5,max(fit$data$data_genotype_non_ref, na.rm=T)), 
     log= 'xy')
abline(b=1, a=0, col = 'black', lty = 2)
axis(1, las = 1, cex.axis = cexaxis, at = c(1, 10, 50, 500), labels = c(1, 10, 50, 500))
axis(2, las = 2, cex.axis = cexaxis, at = c(1, 10, 50, 500), labels = c(1, 10, 50, 500))
for(i in 2:nb_countries){
  points(numbers_obs[,,rd[i]], numbers_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}
cor.test(numbers_obs, numbers_pred_m) ## weight by the number of data point in each pred
## Freqs absolute
freqs_obs = array(NA, dim = c(nb_genotypes, nb_years, nb_countries))
for(i in 1:nb_countries){
  for(j in 1:(nb_genotypes-1)){
    freqs_obs[j,,i] = fit$data$data_genotype_non_ref[j,,i]/fit$data$data_total_number[,i]
  }
  freqs_obs[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]/fit$data$data_total_number[,i]
  
  freqs_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
freqs_obs[which(is.nan(freqs_obs))] = NA
# freqs_obs[which(freqs_obs == 0)] = NA
freqs_pred_m = array(0, dim = c(nb_genotypes, nb_years, nb_countries))
freqs_pred_cimax = array(0, dim = c(nb_genotypes, nb_years, nb_countries))
freqs_pred_cimin = array(0, dim = c(nb_genotypes, nb_years, nb_countries))
for(i in 1:nb_countries){
  count = 1
  t = fit$data$data_total_number[,i]
  for(j in 1:(nb_genotypes-1)){
    tmp = t(apply(Chains$pred_number_non_ref[,i,count,], MARGIN = 1, function(x)x/t))
    freqs_pred_m[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[1])
    freqs_pred_cimin[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[2])
    freqs_pred_cimax[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[3])
    count = count+1
  }
  f = t(apply(Chains$pred_number_ref[,i,], MARGIN = 1, function(x)x/t))
  freqs_pred_m[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[1])
  freqs_pred_cimin[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[2])
  freqs_pred_cimax[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[3])
}
plot(freqs_obs[,,rd[1]], freqs_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), 
     pch = 20, xlab = 'Observed', ylab = 'Predicted', 
     yaxt = 'n', xaxt = 'n', bty='n',
     main = "Frequency",
     cex.axis = cexaxis, cex.lab = cexlab, xlim = c(0,1), ylim = c(0,1))
abline(b=1, a=0, col = 'black', lty = 2)
axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
axis(1, las = 1, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
for(i in 2:nb_countries){
  points(freqs_obs[,,rd[i]], freqs_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}

cor.test(freqs_obs, freqs_pred_m)
# mod = lm(freqs_obs~freqs_pred_m,  weights = numbers_obs_all_clades)
# abline(a = mod$coefficients[1], b = mod$coefficients[2])
# summary(mod)

## Changes predicted numbers
change_obs = fit$data$data_genotype_non_ref[,-1,]/fit$data$data_genotype_non_ref[,-nb_years,]
change_obs[which(change_obs == 0)] = NA
change_pred_m = change_pred_cimin = change_pred_cimax =array(0, dim = c(nb_genotypes-1, nb_years-1, nb_countries))
for(i in 1:nb_countries){
  for(j in 1:(nb_genotypes-1)){
    for(t in 1:(nb_years-1)){
      # change_pred_m[j,t,i] = mean.and.ci(Chains$pred_number[,i,j,t+1]/fit$data$clade_number_paired[j,t,i])[1]
      # change_pred_cimin[j,t,i] = mean.and.ci(Chains$pred_number[,i,j,t+1]/fit$data$clade_number_paired[j,t,i])[2]
      # change_pred_cimax[j,t,i] = mean.and.ci(Chains$pred_number[,i,j,t+1]/fit$data$clade_number_paired[j,t,i])[3]
      change_pred_m[j,t,i] = mean.and.ci(Chains$pred_number_non_ref[,i,j,t+1]/Chains$pred_number_non_ref[,i,j,t])[1]
      change_pred_cimin[j,t,i] = mean.and.ci(Chains$pred_number_non_ref[,i,j,t+1]/Chains$pred_number_non_ref[,i,j,t])[2]
      change_pred_cimax[j,t,i] = mean.and.ci(Chains$pred_number_non_ref[,i,j,t+1]/Chains$pred_number_non_ref[,i,j,t])[3]
    }
  }
  change_obs[,which(colSums(numbers_obs_all_clades[,-1,i]) < threshold),i] = NA
}
plot(change_obs[,,rd[1]], change_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), pch = 20, xlab = 'Observed change', ylab = 'Predicted change', yaxt = 'n', bty = 'n',
     main = "Change (Wrightian fitness)",
     cex.axis = cexaxis, cex.lab = cexlab, xlim = c(0.1,20), ylim = c(0.1,20), cex = cexdots, log = 'xy')
# arrows(change_obs[,,1], change_pred_cimin[,,1],  change_obs[,,1], change_pred_cimax[,,1], length=0.05, angle=90, code=3, 
#        col = 'black')
abline(b=1, a=0, col = 'black', lty = 2, cex.axis = 1.2)
axis(2, las = 2, cex.axis = cexaxis)
for(i in 2:nb_countries){
  points(change_obs[,,rd[i]], change_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}


tmp = change_pred_m
change_pred_m[which(is.infinite(change_pred_m))] = NA
change_obs[which(is.infinite(change_obs))] = NA
cor.test(change_obs, change_pred_m)

mod = lm(change_obs~change_pred_m,  weights = numbers_obs_all_clades[-nb_genotypes,-1,])
summary(mod)

dev.off()
#################################################################################################################################



