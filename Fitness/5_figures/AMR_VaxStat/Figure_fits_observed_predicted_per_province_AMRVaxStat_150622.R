## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(data.table)
## Load fit
setwd("./Fitness/")
fit = readRDS(file = './4_run_model/AMR_VaxStat/output/Output_VaxStatAMR10_fit_all.rds')

## Chains
Chains=rstan::extract(fit$fit)

nb_countries = fit$data$nb_countries
nb_genotypes = fit$data$nb_genotypes_total
nb_years = fit$data$nb_years
nb_countries = fit$data$nb_countries

# ref_clade = 2
numbers_obs_all_clades = array(NA, dim = c(nb_genotypes, nb_years, nb_countries)) 
for(i in 1:nb_countries){
  numbers_obs_all_clades[1:(nb_genotypes),,i] = fit$data$data_genotype_total[,,i]
  # numbers_obs_all_clades[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]
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


load('./3_processed_data/AMR_VaxStat/AMR_VTNVT/clade_sero_list.RData')


pdf(width = 20/2.54, height = 12/2.54, file = "./5_figures/AMR_VaxStat/Figure_fits_shared_param_AMRVAXstat_overall_perprov.pdf", onefile = T)
par(mfcol = c(4,9), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))

colors = c('chartreuse4', 'firebrick', 'royalblue','orange')

for(c in 1:fit$data$nb_countries){
  titles = c('NVT_Susceptible',"NVT_Resistant","VT_Susceptible",'VT_Resistant')
  count =1
  d = fit$data$data_genotype_total[count,,c]
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    count=i
    t = fit$data$data_total_number_total[,c]
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    # if(i < nb_genotypes){
      d = fit$data$data_genotype_total[count,,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
      # count = count +1
    # }
    zeros = fit$data$non_zero_country_year_total[,c]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = Chains$pred_absolute_freq_total[,c,i,1:nb_years]
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    # if(i == 1 & c < 18)  { 
      print(' yay')
      plot(0:(nb_years-1), d_m, type="p", pch=c(rep(15, pch_times-1), rep(15, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,15),cex = cexdots,
           main = paste0(province[c], '-',titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
      # }
    
    
    arrows(0:(nb_years-1),d_ci[1,], 0:(nb_years-1),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    # axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
    axis(1, at = c(1,5, 10, 15), c(0,5, 10, 15)+min_date, cex.axis = cexaxis, lwd = 0.9)
    
    lines(0:(nb_years-1),f_mean_ci[1,], lwd = 1.25,
          cex=1.4,xlab="",ylab="", 
          col = adjustcolor(colors[i], alpha.f = 0.7))
    polygon(x = c(0:(nb_years-1), rev(0:(nb_years-1))),
            y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colors[i], alpha.f = 0.4), border = F)
    
    # abline(v = fit$data$yearF0[c], lty = 3)
    abline(v = fit$data$yearIntroduction[c]-1, lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
}
dev.off()
############################################################################################

############################################################################################
## Plot fits for model  OVERALL
############################################################################################
## Load fit

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6

province = c("Eastern Cape", "Free State", "Gauteng", "KwaZulu-Natal", "Limpopo", "Mpumalanga", "North West", "Northern Cape", "Western Cape")

load('./3_processed_data/AMR_VaxStat/AMR_VTNVT/clade_sero_list.RData')
pdf(width = 12, height = 3, file = "./5_figures/AMR_VaxStat/Figure_fits_shared_param_AMRVAXstat_overall.pdf", onefile = T)
par(mfcol = c(1,4), oma = c(1,1,1,0), mai = c(0.2,0.2,0.2,0.1))
nb_genotypes_total<-fit$data$nb_genotypes_total
colors = c('chartreuse4', 'firebrick', 'royalblue','orange')
fit$data$vaccine_introduction<-rep(10,9)
for(c in 1:1){
  titles = c('NVT_Susceptible',"NVT_Resistant","VT_Susceptible",'VT_Resistant')
  count = 1
  d = rowSums(fit$data$data_genotype_total[count,,])
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  fvt.mat<-d_m.mat<-matrix(nrow=nb_genotypes_total,ncol=nb_years)
  mat.list<-list()
  
  for(i in 1:(nb_genotypes_total)){
    
    mat<-matrix(nrow=nb_years,ncol=7)
    t = rowSums(fit$data$data_total_number_total)
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    # if(i == nb_genotypes_total){
      d = rowSums(fit$data$data_genotype_total[i,,])
      # d = apply(fit$data$data_genotype_total,2,sum)
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
          # tmp = binom.confint(x = d[j], n = t[j], method = c("exact"), type="highest")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
    # }
    zeros = fit$data$non_zero_country_year[,1]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    ##### save the proportion of each in each year
    d_m.mat[i,]<-d_m
    
    f = matrix(0, nrow = nb_chains, ncol = length(1:nb_years))
      for(ccc in 1:nb_countries){
        f = f + t(apply(Chains$pred_absolute_freq_total[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*sum(fit$data$data_total_number_total[,ccc])/sum(t)))
      }

    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    ###store the fits across years for each type
    fvt.mat[i,]<-f_mean_ci[1,]
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    
    
    mat[,1]<-1999+1:nb_years
    mat[,2]<-t(d_m)
    mat[,3:4]<-t(d_ci)
    mat[,5:7]<-t(f_mean_ci)
    colnames(mat)<-c("years","data","datalower","dataupper","fit","fitlower","fitupper")
    mat<-data.table(mat)
    mat$type<-titles[i]
    
    
    plot(1:nb_years, d_m, type="p", pch=c(rep(15, pch_times-1), rep(15, nb_years-pch_times+1)), bty = 'n',
         xlab="Time", ylab = '   ', ylim = ylims,xlim = c(0,15),cex = cexdots,
         main = titles[i], cex.main = cexmain,
         col = adjustcolor('grey30', alpha.f = 0.6), 
         yaxt = 'n', xaxt = 'n', cex.lab = cexlab)
    
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
    abline(v = fit$data$yearIntroduction[c]-1, lty = 3)
    
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
    mat.list[[i]]<-mat
    
  }
}
dev.off()
save(mat.list,file="./5_figures/AMR_VaxStat/overall.fits.RData")
############################################################################################
# Chains$pred_number_total[1,1,1,]
# fit$data$data_genotype_total[1,,1]
# plot(Chains$pred_number_total[1,1,1,],fit$data$data_genotype_total[1,,1])
# plot(Chains$pred_number_total[1,9,2,],fit$data$data_genotype_total[2,,9])
# plot(Chains$pred_number_total[1,1,3,],fit$data$data_genotype_total[3,,1])
# plot(Chains$pred_number_total[1,1,4,],fit$data$data_genotype_total[4,,1])
# 
# abline(b=1,a=0)

