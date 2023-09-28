###fits from NVT and VT parallel AMR plots for Penicillin
## Predicted number vs. observed plot
## Library
library(rstan)
library(RColorBrewer)
library(binom)
library(data.table)

setwd("./Fitness/")
## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

############################################################################################
## Plot options to plot all AMR models with a swithc, without a switch, and without a growth rate
############################################################################################

cexlab = 1.1
cexaxis = 1
cexdots= 1
cexmain = 1.2

province = c("Western Cape","Gauteng","KwaZulu-Natal","Eastern Cape", "Free State","Northern Cape", "North West","Limpopo","Mpumalanga")


###########################################################################################
## LOAD FITS 
############################################################################################
## Load fit
fit = readRDS(file = './4_run_model/AMR_VaxStat/output/Output_vaxstatAMR_NVT.R_fit_all.rds')

## Chains
Chains=rstan::extract(fit$fit)

nb_countries = fit$data$nb_countries
nb_genotypes = fit$data$nb_genotypes_amr
nb_years = fit$data$nb_years

ref_clade = 2

threshold = 5

min_date = 2000
nb_chains = length(Chains$lp__)
############################################################################################

############################################################################################
## Plot fits for model with PENR NVT OVERALL 2011 switch#####
############################################################################################

pdf(width = 3*2, height = 3, file = "./5_figures/AMR_VaxStat/Figure_fits_AMR_NVT_2010_overall.pdf", onefile = T)
par(mfcol = c(1,2), oma = c(1,1,1,0), mai = c(0.2,0.4,0.2,0.4))

colors = c('chartreuse4', 'firebrick')
mat.list<-list()
for(c in 1:1){
  # titles = c("nvt_s" ,"nvt_r", "vt_s", "vt_r")
  titles = c("Susceptible" ,"Resistant")
  
  count = 1
  d =rowSums(fit$data$data_genotype_total[1,,])/rowSums(colSums(fit$data$data_genotype_total[1:2,,]))
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    mat<-matrix(nrow=nb_years,ncol=8)
    t = rowSums(colSums(fit$data$data_genotype_total[1:2,,]))
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    # if(i < nb_genotypes){
      d = rowSums(fit$data$data_genotype_total[i,,])
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
      count = count +1
    # }
    zeros = fit$data$non_zero_country_year_amr[,1]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = matrix(0, nrow = nb_chains, ncol = length(1:nb_years))
      for(ccc in 1:nb_countries){

        tot_predFreq_nvt<-Chains$pred_absolute_freq_total[,ccc,1,1:nb_years]+Chains$pred_absolute_freq_total[,ccc,2,1:nb_years]
        f = f + t(apply(Chains$pred_absolute_freq_total[,ccc,i,1:nb_years]/tot_predFreq_nvt, MARGIN = 1, function(x)x*sum(colSums(fit$data$data_genotype_total[1:2,,])[,ccc])/sum(t)))
      }
  
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print('yay')
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
    
    mat[,1]<-1:nb_years
    mat[,2]<-d_m
    mat[,3:4]<-t(d_ci)
    mat[,5:7]<-t(f_mean_ci)
    mat[,8]<-if(i==1){"Susceptible"}else{"Resistant"}
    mat<-data.table(mat)
    colnames(mat)<-c("years","data","di_lower","di_upper","fit","fit_lower","fit_upper","amr")
    mat.list[[i]]<-mat
    # abline(v = fit$data$yearF0[c], lty = 3)
    # abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
  
}
dev.off()
mat.nvt<-rbindlist(mat.list)
############################################################################################

############################################################################################
## Plot fits for model with PENR VT overall 2011 switch####
############################################################################################

# ## Chains
pdf(width = 3*2, height = 3, file = "./5_figures/AMR_VaxStat/Figure_fits_AMR_VT_2010_overall.pdf", onefile = T)
par(mfcol = c(1,2), oma = c(1,1,1,0), mai = c(0.2,0.4,0.2,0.4))

colors = c('chartreuse4', 'firebrick')
mat.list<-list()
for(c in 1:1){
  # titles = c("nvt_s" ,"nvt_r", "vt_s", "vt_r")
  titles = c("Susceptible" ,"Resistant")
  
  count = 1
  d =rowSums(fit$data$data_genotype_total[1,,])/rowSums(colSums(fit$data$data_genotype_total[3:4,,]))
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    mat<-matrix(nrow=nb_years,ncol=8)
    t = rowSums(colSums(fit$data$data_genotype_total[3:4,,]))
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    # if(i < nb_genotypes){
    d = rowSums(fit$data$data_genotype_total[i+2,,])
    for(j in 1:nb_years){
      if(t[j]>0) {
        tmp = binom.confint(x = d[j], n = t[j], method = c("bayes"), type="central")
        d_m[j] = tmp$mean
        d_ci[1,j] = tmp$lower
        d_ci[2,j] = tmp$upper
      }
    }
    count = count +1
    # }
    zeros = fit$data$non_zero_country_year_amr[,1]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = matrix(0, nrow = nb_chains, ncol = length(1:nb_years))
    # if(i < nb_genotypes){
    for(ccc in 1:nb_countries){
      # f = f + Chains$pred_number_non_ref[,ccc,i,1:nb_years]
      # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*fit$data$data_total_number[,ccc]/t))
      # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*1/9)) 
      tot_predFreq_vt<-Chains$pred_absolute_freq_total[,ccc,3,1:nb_years]+Chains$pred_absolute_freq_total[,ccc,4,1:nb_years]
      f = f + t(apply(Chains$pred_absolute_freq_total[,ccc,i+2,1:nb_years]/tot_predFreq_vt, MARGIN = 1, function(x)x*sum(colSums(fit$data$data_genotype_total[3:4,,])[,ccc])/sum(t)))
    }
    # }
    # if(i == nb_genotypes){
    #   for(ccc in 1:nb_countries){
    #     # f = f + Chains$pred_number_ref[,ccc,1:nb_years]
    #     # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*fit$data$data_total_number[,ccc]/t))
    #     # f = f + t(apply(Chains$pred_absolute_freq[,ccc,i,1:nb_years], MARGIN = 1, function(x)x*1/9))
    #     f = f + t(apply(Chains$pred_absolute_freq_amr[,ccc,i,1:nb_years,1], MARGIN = 1, function(x)x*sum(fit$data$data_total_number_amr[,ccc])/sum(t)))
    #   }
    #   # f = t(apply(f, MARGIN = 1, function(x)x/t))
    # }
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$vaccine_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print('yay')
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
    
    mat[,1]<-1:nb_years
    mat[,2]<-d_m
    mat[,3:4]<-t(d_ci)
    mat[,5:7]<-t(f_mean_ci)
    mat[,8]<-if(i==1){"Susceptible"}else{"Resistant"}
    mat<-data.table(mat)
    colnames(mat)<-c("years","data","di_lower","di_upper","fit","fit_lower","fit_upper","amr")
    mat.list[[i]]<-mat
    # abline(v = fit$data$yearF0[c], lty = 3)
    # abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
  
}
dev.off()
mat.pcv<-rbindlist(mat.list)

mat.pcv$vaxstat<-"VT"
mat.nvt$vaxstat<-"NVT"
df.vaxstatamr<-rbind(mat.pcv,mat.nvt)
setwd("./4_run_model/AMR_VaxStat/")
save(df.vaxstatamr,file="df.vaxstatamr.newmodel.RData")
############################################################################################








