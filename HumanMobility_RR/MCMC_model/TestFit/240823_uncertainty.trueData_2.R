library(dplyr)
library(data.table)
library(ape)
###Incorporating tree uncertainty into data
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/GPS/GPS_GPSC_overall.SA.RData")
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/GPS/GPS_GPSC_everything.RData")

load("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/data/gps_metadata/GPS_GPSC_overall.SA.RData")
load("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/data/gps_metadata/GPS_GPSC_everything.RData")
# ## #Read in population data for each province
load("/Users/sb62/Documents/Migration/SA_Migration_110422/RData_outputs/pop2019_municipality.2017LS.RData")
load("/Users/sb62/Documents/Migration/SA_Migration_110422/RData_outputs/cdr.mat.one.RData")
load('/Users/sb62/Documents/Migration/SA_Migration_110422/RData_outputs/pairwise_geodist.RData')

cdr.mat<-cdr.mat.one
nloc.munic=234
provs<-provs<-colnames(cdr.mat)
splitnames <- matrix(nrow=nloc.munic,ncol=2)
for (us in 1:nloc.munic) {
  for (prov in 1:2) { 
    splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
  }
}

pop_2019<-vector(mode="numeric",length=9)
for (sp in 1:length(provs)){
  sel_sp <- which(splitnames[,2]==provs[sp])
  pop_2019[sp]<-sum(pop2019.town[sel_sp])
}
names(pop_2019)<-provs
## Create dat.in dataframe
nod.quants<-seq(0,1,0.0275)
# nod.quants<-c(0.025,0.5,0.975)

# nboot=length(nod.quants)
### estimating mean and standard deviation from quantiles
# norm_from_quantiles = function(lower, upper, p = 0.25) {
#   mu = mean(c(lower, upper))
#   sigma = (lower - mu) / qnorm(p)
#   list(mu = mu, sigma = sigma)
# }
nboot=50
cummeanDist.emp.bs<-meanDist.emp.bs<-matrix(NaN,101,nboot)
for (j in 1:nboot){
gubgpscs <- c(1,2,5,10,13,14,17,68,79)
dat.GPSC <- output.totEvoTime.SA<-output.typesSA <- outputTipDatesSA <- TimeMRCA.SA <- output.cellindexSA <- output.meanGens.SA <- list()
for (s in gubgpscs) { 
  ##read in all GPSCs
  resbd <- readRDS(paste0("/Users/sb62/Documents/Migration/BactDating/BactDating_Trees_R/bd_GPSC",s))
  ##########sample from the confidence intervals of the nodes##################################################
  # nod.times<-resbd$CI
  # nod.dates<-((resbd$CI[,2]-resbd$CI[,1])/2)+resbd$CI[,1]
  ####sample node date from runif and calculate time to tips
  # node.date<-rep(NA,nrow(nod.times))
  # for(k in 1:nrow(nod.times)){
  #   node.date[k]<-runif(1,min=nod.times[k,][1],max=nod.times[k,][2])
  # }
  

  

  ####sample node date from across CIs and calculate time to tips
  # node.date<-coef_var<-rep(NA,nrow(nod.times))
  # for(k in 1:nrow(nod.times)){
  #   # node.date[k]<-quantile(nod.times[k,],probs=nod.quants[j])
  #   # node.date[k]<-runif(1,min=nod.times[k,][1],max=nod.times[k,][2])
  #   mu<-mean(c(nod.times[k,][1],nod.times[k,][2]))
  #   sigma<-(nod.times[k,][2]-nod.times[k,][1])/4
  #   node.date[k]<-rnorm(1,mu,sigma)
  #   # coef_var[k]<-(norms$sigma/norms$mu)*100
  # }
  # # coef_var_gpscs[[]]<-coef_var
  # nod.dist<-abs(outer(node.date,node.date,"-"))
  # new.edge.length<-rep(NA,length(node.date))
  # edges<-resbd$tree$edge
  # new.edge.length<-mapply(function(i,j) nod.dist[i,j], i=edges[,1],j=edges[,2])
  # ####need edge, Nnode, tip.label
  # samp.tree<-resbd$tree
  # samp.tree$edge.length<-new.edge.length
  
  
  #######Using posterior trees
  Ntips<-Ntip(resbd$tree)
  posts<-resbd$record
  Nposts<-nrow(posts)
  tree_ID = sample(1:Nposts,1) ## I just chose a random number
  tip_dates = resbd$record[tree_ID,1:Ntips]
  node_dates = resbd$record[tree_ID,(Ntips+1):(2*Ntips-1)]
  tipandnodes_dates = c(tip_dates, node_dates)
  node.date<-tipandnodes_dates
  reconstructed_tree = resbd$tree ## Take BD tree for template
  ## Replace edge lengths, based on the tip and node dates that were extracted above
  new.edge.length<-tipandnodes_dates[reconstructed_tree$edge[,2]] - tipandnodes_dates[reconstructed_tree$edge[,1]]
  reconstructed_tree$edge.length = tipandnodes_dates[reconstructed_tree$edge[,2]] - tipandnodes_dates[reconstructed_tree$edge[,1]]
  
  res17<-reconstructed_tree
  ############################################################################################################
  # res17 <- resbd$tree
  lanes <- res17$tip.label
  lanes_tree <- subset(lanes, lanes  %in% GPS_GPSC_overall$Lane_Id)
  trees <- keep.tip(res17, lanes_tree)
  #Reorder tree based on order of tip labels
  GPS_GPSC17 <- subset(GPS_GPSC_overall, GPS_GPSC_overall$Lane_Id %in% lanes)
  order <- trees$tip.label 
  GPS_GPSC17 <- GPS_GPSC17 %>%
    dplyr::slice(match(order, Lane_Id))
  #TipIDs
  labels<- trees$tip.label
  labels <- 1:length(labels)
  label.dat <- outer(labels,labels,"==")
  rownames(label.dat) <- labels
  colnames(label.dat) <- labels
  tipIDS <- data.frame(tip1 = as.integer(),tip2 = as.integer())
  tipIDS <- suppressWarnings(melt(label.dat)[,c(1:2)])
  #Locations
  provs<-GPS_GPSC17$Region
  tmp = provs
  tmp[which(tmp=="Eastern Cape")] <- 1
  tmp[which(tmp=="Free State")] <- 2
  tmp[which(tmp=="Gauteng")] <- 3
  tmp[which(tmp=="KwaZulu-Natal")] <- 4
  tmp[which(tmp=="Limpopo")] <- 5
  tmp[which(tmp=="Mpumalanga")] <- 6
  tmp[which(tmp=="North West")] <- 7
  tmp[which(tmp=="Northern Cape")] <- 8
  tmp[which(tmp=="Western Cape")] <- 9
  provs.dat <- outer(tmp,tmp,"==")
  rownames(provs.dat) <- provs
  colnames(provs.dat) <- provs
  locations <- data.frame(loc1=as.numeric(),loc2=as.numeric())
  locations <- suppressWarnings(melt(provs.dat)[,c(1:2)])
  #MRCAs
  dates<-GPS_GPSC17$Col_time
  dist.mat<-cophenetic.phylo(trees)
  pat.dist.matGPSC<-matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))
  diag(pat.dist.matGPSC)<-NA
  mrca<-(outer(dates,dates,"+")-pat.dist.matGPSC)/2
  mrcas <- data.frame(mrca=as.numeric())
  mrcas <- suppressWarnings(melt(mrca)[,c(3)])
  totEvoTime <- suppressWarnings(melt(pat.dist.matGPSC)[,c(3)])
  
  
  dat.GPSC[[which(gubgpscs==s)]]<-cbind(tipIDS,locations,mrcas,totEvoTime)
  colnames(dat.GPSC[[which(gubgpscs==s)]]) <- c("tip1","tip2","loc1","loc2","mrca","totEvoTime")
  dat.GPSC[[which(gubgpscs==s)]]$tip1 <- as.integer( dat.GPSC[[which(gubgpscs==s)]]$tip1)
  dat.GPSC[[which(gubgpscs==s)]]$tip2 <- as.integer( dat.GPSC[[which(gubgpscs==s)]]$tip2)
  dat.GPSC[[which(gubgpscs==s)]]$loc1 <- as.numeric( dat.GPSC[[which(gubgpscs==s)]]$loc1)
  dat.GPSC[[which(gubgpscs==s)]]$loc2 <- as.numeric( dat.GPSC[[which(gubgpscs==s)]]$loc2)
  dat.GPSC[[which(gubgpscs==s)]]$mrca <- as.numeric( dat.GPSC[[which(gubgpscs==s)]]$mrca)
  dat.GPSC[[which(gubgpscs==s)]]$totEvoTime <- as.numeric( dat.GPSC[[which(gubgpscs==s)]]$totEvoTime)
  
  ### Carriage or Disease
  types <- GPS_GPSC17$Type
  output.typesSA[[which(gubgpscs==s)]] <- types
  #Additional necessary dataframes
  outputTipDatesSA[[which(gubgpscs==s)]] <- dates #Sample Isolation dates 
  TimeMRCA.SA[[which(gubgpscs==s)]] <- -sweep(mrca,1,dates,"-")
  output.totEvoTime.SA[[which(gubgpscs==s)]] <- pat.dist.matGPSC*365
  
  tmp = provs
  tmp[which(tmp=="Eastern Cape")] <- 1
  tmp[which(tmp=="Free State")] <- 2
  tmp[which(tmp=="Gauteng")] <- 3
  tmp[which(tmp=="KwaZulu-Natal")] <- 4
  tmp[which(tmp=="Limpopo")] <- 5
  tmp[which(tmp=="Mpumalanga")] <- 6
  tmp[which(tmp=="North West")] <- 7
  tmp[which(tmp=="Northern Cape")] <- 8
  tmp[which(tmp=="Western Cape")] <- 9
  output.cellindexSA[[which(gubgpscs==s)]] <- as.numeric(tmp)
}
names(dat.GPSC) <- gubgpscs 
## List with parameter values
extTranMatDat.GPSC <- list(pars=c(),popbyCell=pop_2019)
dat.tmp.allser<-list()
for (ii in 1:9){
  dat.tmp.allser[[ii]]<-cbind(expand.grid(1:length(output.cellindexSA[[ii]]),1:length(output.cellindexSA[[ii]])),
                              expand.grid(output.cellindexSA[[ii]],output.cellindexSA[[ii]]),
                              as.numeric(output.totEvoTime.SA[[ii]]),
                              expand.grid(outputTipDatesSA[[ii]],outputTipDatesSA[[ii]]))
  colnames(dat.tmp.allser[[ii]])<-c("id1","id2","loc1","loc2","totTimeDays","cy1","cy2")
  a<-which(dat.tmp.allser[[ii]][,1]>dat.tmp.allser[[ii]][,2]&is.na(rowSums(dat.tmp.allser[[ii]][,1:4]))==F)
  
  dat.tmp.allser[[ii]]<-dat.tmp.allser[[ii]][a,]
}




#############################################
### Prepare data
nByGPSC<-rep(NaN,9)
for (i in 1:9)nByGPSC[i]<-nrow(dat.tmp.allser[[i]])
dat.inMaster<-do.call("rbind",dat.tmp.allser)
dat.inMaster<-dat.inMaster[which(is.na(dat.inMaster$totTimeDays)==F),]
dat.inMaster$GPSC<-rep(1:9,nByGPSC)

dat.inMaster$year1<-floor(dat.inMaster$cy1)
dat.inMaster$year2<-floor(dat.inMaster$cy2)
a<-as.matrix(dat.inMaster[,3:4])
dat.inMaster$SDist<-pairwise_geodist[a]


### Sampling probability by year, GPSC and province
NoDetectByYearGPSC<-table(c(dat.inMaster$loc1,dat.inMaster$loc2),c(dat.inMaster$GPSC,dat.inMaster$GPSC),c(dat.inMaster$year1,dat.inMaster$year2))
un.years<-sort(unique(dat.inMaster$year1))

avNoDetectByloc<-table(c(dat.inMaster$loc1,dat.inMaster$loc2))
dat.inMaster$YearREF1<-match(dat.inMaster$year1,un.years)
dat.inMaster$YearREF2<-match(dat.inMaster$year2,un.years)

dat.inMaster$SDist[which(dat.inMaster$SDist==0)]<-10
pairwise_geodist[which(pairwise_geodist==0)]<-10

id.1<-paste0(dat.inMaster$id1,".",dat.inMaster$GPSC)
id.2<-paste0(dat.inMaster$id2,".",dat.inMaster$GPSC)
un.ids<-unique(c(id.1,id.2))
nids<-length(un.ids)

gDistMax<-seq(0,100,1)
gDistMin<-gDistMax-1
gDistMin[which(gDistMin<0)]<-0
gDistMid<-(gDistMin+gDistMax)/2

###Add evolutionary time in years column
dat.inMaster$totEvoTime<-dat.inMaster$totTimeDays/365

for (i in 1:length(gDistMax)){
  a<-which(dat.inMaster$totEvoTime>=gDistMin[i]&dat.inMaster$totEvoTime<gDistMax[i])
  b<-which(dat.inMaster$totEvoTime<gDistMax[i])
  meanDist.emp.bs[i,j]<-mean(dat.inMaster$SDist[a])
  cummeanDist.emp.bs[i,j]<-mean(dat.inMaster$SDist[a])
}
print(paste0("Bootstrap: ",j))
}
cummeanDist.emp.bs.qt<-apply(cummeanDist.emp.bs,1,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
meanDist.emp.bs.qt<-apply(meanDist.emp.bs,1,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
overall_data<-matrix(nrow=length(gDistMid),ncol=4)
overall_data[,1]<-gDistMid
overall_data[,2:4]<-t(meanDist.emp.bs.qt)
overall_data<-data.table(overall_data)
colnames(overall_data)<-c("years","lowerCI","median","upperCI")

ggplot()+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=overall_data,aes(x=years,y=median),color="blue")+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),
        legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
  geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="blue")+
  xlab("Evolution Time (years)")+
  ylab("Distance (km)")#+
  xlim(0,50)

# setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/TestFit/")
save(overall_data,file="./Descriptive_Stats_plots/overall_data_multiTree.RData")
save(overall_data,file="./MCMC_model/TestFit/Data/overall_data_multiTree.RData")



  
  
  #########################################################
  ####Calculate the coefficient of variation for each node  ###  ###  ###  ###
#   ############################################################
#   gubgpscs <- c(1,2,5,10,13,14,17,68,79)
#   coef_var_gpscs<-coef_var_gpscs_lt10<-list()
#   for (s in gubgpscs) { 
#     ##read in all GPSCs
#     resbd <- readRDS(paste0("/Users/sb62/Documents/Migration/BactDating/BactDating_Trees_R/bd_GPSC",s))
#     ##########sample from the confidence intervals of the nodes##################################################
#     nod.times<-resbd$CI
#     ### find median node date from the CIs
#     nod.dates<-((resbd$CI[,2]-resbd$CI[,1])/2)+resbd$CI[,1]
#     ### identify where the nodes that aren't leaves start
#     tipnum<-length(resbd$tree$tip.label)+1
#     rownums<-tipnum:nrow(nod.times)
#     
#     ### identify all pairs that are <10 years divergent
#     nod.dates_nts<-nod.dates[rownums]
#     nod.diffs<-abs(outer(nod.dates_nts,nod.dates_nts,"-"))
#     samp_lt10<-which(nod.diffs<10)
#     ####calculate the coefficient of variation for all tree tips
#     node.date<-coef_var<-rep(NA,nrow(nod.times)-tipnum)
#     for(k in 1:length(rownums)){
#       mu<-mean(c(nod.times[k,][1],nod.times[k,][2]))
#       sigma<-(nod.times[k,][2]-nod.times[k,][1])/4
#       node.date[k]<-rnorm(1,mu,sigma)
#       ## calculate the coefficient of variation
#       coef_var[k]<-(sigma/mu)*100
#     }
#     print(s)
#     tmp<-outer(coef_var,coef_var,FUN="paste")
#     coef_varlt10<-tmp[samp_lt10]
#     coef_var_gpscs[[which(gubgpscs==s)]]<-coef_var
#     coef_var_gpscs_lt10[[which(gubgpscs==s)]]<-coef_varlt10
#   }
#   
#   
# coef_var_gpscs_num<-list()
#   for (i in 1:9){
#     numpairs<-length(coef_var_gpscs_lt10[[i]])
#     tmp_mat<-matrix(nrow=numpairs,ncol=2)
#     for ( j in 1:numpairs){
#       tmp_mat[j,1]<-as.numeric(str_split(coef_var_gpscs_lt10[[i]][j]," ")[[1]][1])
#       tmp_mat[j,2]<-as.numeric(str_split(coef_var_gpscs_lt10[[i]][j]," ")[[1]][2])
#     }
#     coef_var_gpscs_num[[i]]<-tmp_mat
#     print(i)
#   }
# 
# ###coefficient of variance around all nodes with <10 years divergence to branches
# mat_sum<-matrix(nrow=9,ncol=5)
# for(i in 1:9){
#   vec_coefs<-c(coef_var_gpscs_num[[i]][,1],coef_var_gpscs_num[[i]][,2])
#   mat_sum[i,1:3]<-quantile(vec_coefs,probs=c(0.025,0.5,0.975))
#   mat_sum[i,4]<-mean(vec_coefs)
#   mat_sum[i,5]<-gubgpscs[i]
# }
#  
# colnames(mat_sum)<-c("lowerCI","median","upperCI","mean","GPSC") 
# 
# 
# 
# #########################################################
# ####Calculate the coefficient of variation for each tip to node distance in generations  ###  ###  ###  ###
# ############################################################
# library(data.table)
# library(ape)
# library(dplyr)
# #### calc. coeffecient of variance from node to tip
# # sd(date sample - date node)/mean(date sample-date node)
# # norm_from_quantiles = function(lower, upper, p = 0.25) {
# #   mu = mean(c(lower, upper))
# #   sigma = (lower - mu) / qnorm(p)
# #   list(mu = mu, sigma = sigma)
# # }
# gubgpscs <- c(1,2,5,10,13,14,17,68,79)
# gub_cov<-matrix(nrow=length(gubgpscs),ncol = 5)
# for(g in 1:length(gubgpscs)){
#   resbd <- readRDS(paste0("/Users/sb62/Documents/Migration/BactDating/BactDating_Trees_R/bd_GPSC",gubgpscs[g]))
# ### split CIs into leaves and nodes and drop the deepest node
# nod.times<-resbd$CI
# ### pull tip times
# tip.times<-nod.times[which(nod.times[,1]==nod.times[,2]),1]
# ### pull node times intervals
# intnod.times_CI<-nod.times[which(nod.times[,1]!=nod.times[,2]),]
# ### find the mediam across node CIs
# nod.dates<-((intnod.times_CI[,2]-intnod.times_CI[,1])/2)+intnod.times_CI[,1]
# ### distance from nodes to tips
# mat<-abs(outer(nod.dates,tip.times,"-"))
# mat<-c(mat)
# whichfit<-which(mat<10)
# #### sample from the normal distribution around each node and re-estimate the difference in time to the tips
# nboot=100
# mat_diff_boot<-matrix(nrow=length(whichfit),ncol=nboot)
# for(i in 1:nboot){
#   nod.dates_samp<-vector(mode="numeric",length=length(nod.dates))
#   for (k in 1:nrow(intnod.times_CI)){
#     # repeat{
#     # norms<-norm_from_quantiles(intnod.times_CI[k,][1],intnod.times_CI[k,][2])
#     # nod.dates_samp[k]<-rnorm(1,mean=norms$mu,sd=norms$sigma)
#     
#     mu<-mean(c(intnod.times_CI[k,][1],intnod.times_CI[k,][2]))
#     sigma<-(intnod.times_CI[k,][2]-intnod.times_CI[k,][1])/4
#     # nod.dates_samp[k]<-rnorm(1,mu,sigma)
#     repeat{
#       nod.dates_samp[k]<-rnorm(1,mu,sigma)
#       # print(k)
#       if(nod.dates_samp[k]<=intnod.times_CI[k,][2]){
#         break
#       }
#     }
#     if(100000000%%k==0){print(k)}
#   }
#     vec_diff<-abs(outer(nod.dates_samp,tip.times,"-"))
#     # vec_diff<-vec_diff[upper.tri(vec_diff)]
#     vec_diff<-c(vec_diff)
#     vec_diff<-vec_diff[whichfit] ## only include those with medians which are <10 years from a tip
#     vec_diff<-(vec_diff*365)/35 ### convert them to number of generations
#     mat_diff_boot[,i]<-vec_diff  
#  print(i)  
# }
# mean_vec<-apply(mat_diff_boot,1,mean)
# std_vec<-apply(mat_diff_boot,1,sd)
# cov_vec<-(std_vec/mean_vec)
# gub_cov[g,1]<-gubgpscs[g]
# gub_cov[g,2:4]<-quantile(cov_vec[whichfit],probs=c(0.025,0.5,0.975),na.rm=T)
# gub_cov[g,5]<-mean(cov_vec[whichfit],na.rm=T)
# print(g)
# }
# #
# 
# colnames(gub_cov)<-c("gpsc","lowerCI","median","upperCI","mean")
# 
# 
# gub_cov_lt2<-gub_cov
# 
# gub_cov_lt10<-gub_cov
