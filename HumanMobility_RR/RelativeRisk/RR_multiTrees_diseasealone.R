######################
###########FUNCTIONS
######################
bind.trees<-function(trees){
  if(length(trees)==2) return(bind.tree(trees[[1]],trees[[2]]))
  else { 
    # trees<-c(bind.tree(trees[[1]],trees[[2]]),
             # if(length(trees)>2) trees[3:length(trees)])
    # trees<-bind.trees(trees) ## this is the recursive part
    ## #Bind trees
    res1=trees[[1]]
    res2=trees[[2]]
    res5=trees[[3]]
    res14=trees[[4]]
    res17=trees[[5]]
    res68=trees[[6]]
    res79=trees[[7]]
    res10=trees[[8]]
    res13=trees[[9]]
    resA <- bind.tree(res2,res5)
    resB <- bind.tree(res14,res17)
    resC <- bind.tree(res68, res79)
    resD <- bind.tree(res10, res13)
    resE <- bind.tree(resA,resB)
    resF <- bind.tree(resC,resD)
    resG <- bind.tree(resE, resF)
    res_everything <- bind.tree(resG,res1)
    
    return(res_everything)
  }
}

neighbor_cont_bootstrap <- function(x=x,y=x, geo_matrix, colyear_matrix, strain_matrix, gpsc_matrix,denom_matrix){
  geo_mat.tmp = geo_matrix[x,y]
  time_mat.tmp = colyear_matrix[x,y]
  strain_mat.tmp = strain_matrix[x,y]
  denom_mat.tmp = denom_matrix[x,y]
  gpsc_mat.tmp = gpsc_matrix[x,y]
  
  tmp = (geo_mat.tmp) * time_mat.tmp 
  tmp[which(tmp==0)]<-NA
  
  tmp2 = (geo_mat.tmp*(geo_mat.tmp>1000)*( geo_mat.tmp<1600)) * time_mat.tmp
  tmp2[which(tmp2==0)]<-NA
  
  a1=cumsum(hist(tmp*strain_mat.tmp,breaks=c(0,0,minpos,1e10),plot=F)$counts)
  a2=cumsum(hist(tmp*strain_mat.tmp,breaks=c(0,0,maxpos,1e10),plot=F)$counts)
  a=a2-a1
  
  b1=cumsum(hist(tmp*(gpsc_mat.tmp),breaks=c(0,0,minpos,1e10),plot=F)$counts)
  b2=cumsum(hist(tmp*(gpsc_mat.tmp),breaks=c(0,0,maxpos,1e10),plot=F)$counts)
  b=b2-b1
  
  c = rep(sum(tmp2*strain_mat.tmp, na.rm=T),length(b))
  
  d = rep(sum(tmp2*gpsc_mat.tmp, na.rm=T),length(a))
  
  rr = (a/b)/(c/d) 
  if (max(rr, na.rm = T)>1E10){
    rr = ((a+1)/(b+1))/((c+1)/(d+1))
  }
  if (min(rr, na.rm=T)==0){
    rr = ((a+1)/(b+1))/((c+1)/(d+1)) 
  }
  rr = rr[2:(length(rr)-1)]
  return(rr)
}
######################################
#################################
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(geodist)
library(ape)
library(BactDating)

gpscgub <- c(2,79,1,5,10,13,14,17,68)

categorical=FALSE
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
load("./data/gps_metadata/GPS_SA.RData")

##Load in Data
# load("../data/gps_metadata/GPS_GPSC_everything.RData")
'%notin%'<-Negate('%in%')
load("./data/gps_metadata/GPS_GPSC_everything.RData")
GPS_GPSC_everything.Disease<-subset(GPS_GPSC_everything,GPS_GPSC_everything$Type=="Disease")
nboot=50
# boot.out = matrix(NA, length(maxpos), nboot )
boot_list<-list()
if(categorical==TRUE){
  boot_array<-array(NA,dim=c(6,4,nboot))
}else{
boot_array<-array(NA,dim=c(6,9,nboot))}
for (j in 1:(nboot)){
  # #Read in trees
  trees<-list()
  for(g in gpscgub){
    resbd <- readRDS(paste0("./data/phylogenies/bd_GPSC",g))
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
    trees[[which(gpscgub==g)]]<-reconstructed_tree
  }
  res_everything<-bind.trees(trees)
  lanes <- res_everything$tip.label
  
  
  
# resbd <- readRDS("./data/phylogenies/bd_GPSC2")
# res2 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC5")
# res5 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC14")
# res14 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC17")
# res17 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC13")
# res13 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC10")
# res10 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC68")
# res68 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC79")
# res79 <- resbd$tree
# resbd <- readRDS("./data/phylogenies/bd_GPSC1")
# res1 <- resbd$tree
# ## #Bind trees
# resA <- bind.tree(res2,res5)
# resB <- bind.tree(res14,res17)
# resC <- bind.tree(res68, res79)
# resD <- bind.tree(res10, res13)
# resE <- bind.tree(resA,resB)
# resF <- bind.tree(resC,resD)
# resG <- bind.tree(resE, resF)
# res_everything <- bind.tree(resG,res1)
# lanes <- res_everything$tip.label




# ###### Drop tips of tree not included in metadata
lanes_tree <- subset(lanes, lanes  %in% GPS_GPSC_everything.Disease$Lane_Id)
res_overall <- keep.tip(res_everything, lanes_tree)
# #Reorder tree based on order of tip labels
order <- res_overall$tip.label
GPS_GPSC_everything.Disease <- GPS_GPSC_everything.Disease %>%
  dplyr::slice(match(order, Lane_Id))


#### Create Matrices
## #Create Geo Matrix
vec_coords = cbind(as.numeric(GPS_GPSC_everything.Disease$Longitude), as.numeric(GPS_GPSC_everything.Disease$Latitude))
geo_mat_manyCont= geodist(
  vec_coords,
  vec_coords,
  paired=FALSE,
  measure="haversine")/1000
geo_mat_manyCont <- round(geo_mat_manyCont, 2)
diag(geo_mat_manyCont) <- NA
geo_mat_manyCont[geo_mat_manyCont==0] <- 10
colnames(geo_mat_manyCont) <-  GPS_GPSC_everything.Disease$Region
rownames(geo_mat_manyCont) <-  GPS_GPSC_everything.Disease$Region
## #Create similarity matrix
dist.mat <- cophenetic.phylo(res_overall)
diag(dist.mat) <- NA
## Create collection year matrix
lane.names <- GPS_GPSC_everything.Disease$Lane_Id
vector_colyear <- GPS_GPSC_everything.Disease$Col_time #Exact collection years
colyear_mat = abs(outer(vector_colyear, vector_colyear, "-"))
rownames(colyear_mat) <- lane.names
colnames(colyear_mat) <- lane.names
diag(colyear_mat) <- NA
strain_mat <- (dist.mat-colyear_mat)/2

## Find distances to other countries
lat_vec <- unique(as.numeric(GPS_GPSC_everything.Disease$Latitude))
lat_vec <- subset(lat_vec,!is.na(lat_vec))
long_vec <- unique(as.numeric(GPS_GPSC_everything.Disease$Longitude))
long_vec <- subset(long_vec,!is.na(long_vec))
region_vec <- unique(GPS_GPSC_everything.Disease$Region)
x <- cbind(long_vec,lat_vec)
y <- x
pairwise_distCont <- geodist(
  x,
  y,
  paired = FALSE,
  measure = "haversine") /1000
colnames(pairwise_distCont) <- region_vec
rownames(pairwise_distCont) <- region_vec
pairwise_distCont<- round(pairwise_distCont, 2)

##Create vector of pairwise distances for each continent and assign specific pairwise distances for within Africa and beyond Africa

outsideAfrica <- unique(GPS_GPSC_everything.Disease$Country[GPS_GPSC_everything.Disease$Continent != "Africa"])
outsideAfrica_distvec <- pairwise_distCont[rownames(pairwise_distCont) %in% outsideAfrica][pairwise_distCont[rownames(pairwise_distCont) %in% outsideAfrica]!=0]

outsideSouthAfrica <- unique(GPS_GPSC_everything.Disease$Country[GPS_GPSC_everything.Disease$Continent == "Africa" & GPS_GPSC_everything.Disease$Country != "South Africa"])
outsideSouthAfrica_dist_vec <-pairwise_distCont[rownames(pairwise_distCont) %in% outsideSouthAfrica][pairwise_distCont[rownames(pairwise_distCont) %in% outsideSouthAfrica]!=0]

geo_mat_manyCont[geo_mat_manyCont %in% outsideSouthAfrica_dist_vec] <- 2500
geo_mat_manyCont[geo_mat_manyCont %in% outsideAfrica_distvec] <- 6500


##Extra Matrices
same_country <- outer(GPS_GPSC_everything.Disease$Country,GPS_GPSC_everything.Disease$Country,"==")
diag(same_country) <- NA
same_country = same_country+0
## GPSC
gpsc_mat <- outer(GPS_GPSC_everything.Disease$GPSC,GPS_GPSC_everything.Disease$GPSC,"==")
diag(gpsc_mat) <- NA
gpsc_mat <- gpsc_mat + 0
vector_gpsc <- GPS_GPSC_everything.Disease$GPSC
# Weight vector for sampling
weight_mat <-data.table(table(GPS_GPSC_everything.Disease$Region, GPS_GPSC_everything.Disease$GPSC))
colnames(weight_mat) <- c("Region","GPSC","N")
dataset_new <- GPS_GPSC_everything.Disease %>% left_join(weight_mat, by=c("Region","GPSC"))
weight_vector <- dataset_new$N
if(categorical==TRUE){
# #To look at different collection time ranges
tseqmax=c(5,10,20,200)
tseqmin=c(0,5,10,20)}else{
## to look continuously across time ranges
tseqmax=seq(1,90,10)
tseqmin=tseqmax-20
tseqmin[which(tseqmin<0)] <- 0
}
plot_roll <- list()
mat_list <- list()
mat_plot <- list()
mat_country <-  vector(mode = "list", length = 2)
vector_region <- GPS_GPSC_everything.Disease$Region
# dists matrix
maxpos=c(100,500,1000,1600,c(3000,7000))
minpos=c(0,100,500,1000,1990,4000)
tseq_out<-matrix(ncol=length(tseqmax),nrow=length(maxpos))
for (k in 1:length(tseqmax)) {
  sa_row <- which(GPS_GPSC_everything.Disease$Country == "South Africa" )
  strain_matrix = (strain_mat>tseqmin[k])*(strain_mat<=tseqmax[k])
  geo_matrix = geo_mat_manyCont
  colyear_matrix = (colyear_mat>0)*(colyear_mat<=1)
  gpsc_matrix = gpsc_mat
  denom_matrix = (strain_mat>tseqmax[k])
  nseq = dim(colyear_matrix)[1]
  region_vec <- names(table(GPS_GPSC_everything.Disease$Region))
  ### multiple samples across samples from regions
  sample_boot<-vector(mode="numeric",length=20)
  for (sb in 1:20){
      selectedSeqs=NULL
      for (i in 1:length(region_vec)) {
        a <- which(vector_region==region_vec[i] )
        b <- if (length(a) == 1) a else sample(a,min(300,length(a)),replace=F)
        selectedSeqs <- c(selectedSeqs,b)
      }
      # selectedSeqs <- sample(nseq,replace=T,prob=c(1/weight_vector))#length of number of sequences of that location with that GPSCs
      tmp = sample(selectedSeqs, replace = T)
    rr =neighbor_cont_bootstrap(x=sa_row,y=tmp,geo_matrix=geo_matrix, colyear_matrix=colyear_matrix, strain_matrix=strain_matrix, gpsc_matrix=gpsc_matrix,denom_matrix=denom_matrix)
    # boot.out[,j] = rr
    sample_boot[sb]<-rr
  }
    print(k)
    # tseq_out[,k]<-rr
    tseq_out[,k]<-mean(sample_boot)
    
    # print(j)
}

boot_list[[j]]<-tseq_out
boot_array[,,j]<-tseq_out
}

# tseq_out<-data.table(tseq_out)
# colnames(tseq_out)<-paste(tseqmin,"-",tseqmax)
# rownames(tseq_out)<-paste(minpos,"-",maxpos)
ts_list<-list()
for(ts in 1:length(tseqmax)){
  mat<-matrix(nrow=length(maxpos),ncol=7)
  for(mp in 1:length(maxpos)){
    boot_vec<-boot_array[mp,ts,]
    mat[mp,1] <- paste(minpos[mp],"-",maxpos[mp])
    mat[mp,2:4] <- quantile(boot_vec,probs=c(0.5,0.025,0.975))
    # mat[,3] <- boot.ci[1,] 
    # mat[,4] <- boot.ci[2,] 
    mat[mp,5] <- paste(tseqmin[ts],"-",tseqmax[ts])
    mat[mp,6] <- (tseqmax[ts]+tseqmin[ts])/2
    mat[mp,7] <- minpos[mp]
  }
  mat<-data.table(mat)
  ts_list[[ts]]<-mat
}

mat_sum<-rbindlist(ts_list)
  ## Matrix
  colnames(mat_sum) <- c("distance_range" , "RR", "lowerCI", "upperCI", "time_range", "medMRCA","minimum_dist")
  mat_sum <- as.data.table(mat_sum)
  mat_sum$RR <- as.numeric(mat_sum$RR)
  mat_sum$lowerCI <- as.numeric(mat_sum$lowerCI)
  mat_sum$upperCI <- as.numeric(mat_sum$upperCI)
  mat_sum$minimum_dist <- as.numeric(mat_sum$minimum_dist)
  mat_sum$medMRCA <- as.numeric(mat_sum$medMRCA)
  mat_sum$distance_range[mat_sum$distance_range=="0 - 100"]<-"Within Province"
  mat_sum$distance_range[mat_sum$distance_range=="100 - 500"]<-"<500"
  mat_sum$distance_range[mat_sum$distance_range=="500 - 1000"]<-"500-1000"
  mat_sum$distance_range[mat_sum$distance_range=="1000 - 1600"]<-"Distant Pairs"
  mat_sum$distance_range[mat_sum$distance_range=="1990 - 3000"]<-"Other Africa"
  mat_sum$distance_range[mat_sum$distance_range=="4000 - 7000"]<-"Outside Africa"
  
  
  
  
  
  
  
  ####################################
  # PLOTS
if(categorical==TRUE){
  matplot_categorical<-mat_sum
}else{
matplot_continuous <- mat_sum
}

if(categorical==TRUE){
  load("./RelativeRisk/files/mat_lineage.Sub.Rdata")###This loads all 6910 SA metadata
  mat_lineage.noSub<-mat_lineage
  mat_lineage.noSub$time_range[which(mat_lineage.noSub$time_range=="0 - 1")]<-"Same Lineage"
  
  matplot_specificRange.noSub<-matplot_categorical
  save(matplot_specificRange.noSub,file="./RelativeRisk/files/matplot_specificRange.noSub_cat_mulTrees_260923_disease.RData")
  
  
}else{
  save(matplot_continuous,file="./RelativeRisk/files/mat_sub_genomic_cont_multiTrees_disease.RData")
  
}
  
if(categorical==TRUE){
  load("./RelativeRisk/files/matplot_specificRange.noSub_cat_mulTrees_260923_disease.RData")
  mat.tmp <- rbind(matplot_specificRange.noSub[,c("distance_range","RR","lowerCI","upperCI","time_range")],
                   mat_lineage.noSub[,c("distance_range","RR","lowerCI","upperCI","time_range")])
  mat.tmp$lowerCI[mat.tmp$lowerCI <= 0.1 ] <- 0.1
  mat.tmp$upperCI[mat.tmp$upperCI <= 0.1]<- 0.1
  mat.tmp$RR[mat.tmp$RR <= 0.1]<- 0.1
  mat.tmp$lowerCI[mat.tmp$lowerCI >= 10 ] <- 10
  mat.tmp$upperCI[mat.tmp$upperCI >= 10 ] <- 10
  mat.tmp$RR[mat.tmp$RR >= 10 ] <- 10
  
  mat.tmp$time_range_f = factor(mat.tmp$time_range, levels=c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"), labels = c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"))
  
  mat.tmp$time_range_f = factor(mat.tmp$time_range, levels=c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"), labels = c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"))
  mat.tmp$distance_range_f = factor(mat.tmp$distance_range, levels=c("Within Province", "<500","500-1000","Distant Pairs" ,"Other Africa","Outside Africa"))
  # pdf(file="/Users/sb62/Documents/Migration/SA_Migration_110422/Submission/Figures/Figure_panel2/RRoverDistTime.PDF.pdf",width=3.4,height=10)
  mat.tmp$upperCI[which(mat.tmp$upperCI>5)]<-5
  lineplot <- ggplot(data = mat.tmp, aes( x = distance_range_f, y = RR , group = time_range_f)) +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    geom_point(size=1) +
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), alpha = .9, color = "black", width = 0.5 ) +
    theme(axis.text.x = element_text(angle = 90, size=20),
          axis.text.y = element_text(size=18), axis.title.y  = element_text(size=25),
          axis.title.x=element_text(size=18), title =element_text(size=18),
          panel.grid.major= element_blank(), panel.grid.minor= element_blank(),axis.line = element_line(colour = "black")) +
    scale_x_discrete(name ="",
                     
                     limits=c("Within Province", "<500","500-1000","Distant Pairs" ,"Other Africa","Outside Africa"))+
    # theme(aspect.ratio=20/30) +
    
    geom_point(aes(x=4, y=1),shape=25,fill="white", size=1,alpha=1)+
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_y_continuous( trans = "log10", breaks = c(0.09,1,5),labels=c("<0.1","1.0","5.0"),limits=c(0.09,5)) +
    ylab("Risk Ratio")
  
  lin1_mal_cont <- lineplot + geom_rect(fill = '#00539CFF', xmin = 0, xmax = 1.5, ymin =-5, ymax = 7, alpha =0.05) +
    geom_rect(fill = '#ED2B33FF', xmin = 1.5, xmax = 4.5, ymin =-5, ymax = 7, alpha =0.05) +   
    geom_rect(fill = '#97BC62FF', xmin = 4.5, xmax = 5.5, ymin =-5, ymax = 7, alpha =0.05) +
    geom_rect(fill = 'purple', xmin = 5.5, xmax = 6.5, ymin =-5, ymax = 7, alpha =0.05) #+
  pmanycont <- lin1_mal_cont + facet_wrap( ~ time_range_f, nrow =5) +
    theme(strip.text.x=element_blank())
  pmanycont
  
  
}else{
load("./RelativeRisk/files/mat_sub_genomic_cont_multiTrees_disease.RData")
  mat_distancerange <- matplot_continuous
  mat_distancerange$lowerCI[mat_distancerange$lowerCI <= 0.1 ] <- 0.1
  mat_distancerange$upperCI[mat_distancerange$upperCI <= 0.1]<- 0.1
  mat_distancerange$RR[mat_distancerange$RR <= 0.1]<- 0.1
  
  mat_sub_genomic <- subset(mat_distancerange, mat_distancerange$distance_range == "Within Province" | mat_distancerange$distance_range == "500-1000" | mat_distancerange$distance_range == "Other Africa"| mat_distancerange$distance_range == "Outside Africa")
  mat_sub_genomic$distance_range_f = factor(mat_sub_genomic$distance_range, levels= c("Within Province","500-1000","Other Africa","Outside Africa"))
  mat_sub_genomic.tmp <- subset(mat_sub_genomic, mat_sub_genomic$distance_range=="Within Province")
  }
plotwithin <- ggplot(data = mat_sub_genomic, aes( x = medMRCA, y = RR, group = distance_range_f ,color = distance_range_f)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI), alpha = 0.1, color = NA, fill = "#00539CFF" ) +
  theme(axis.text.x = element_text(angle = 90, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        title = element_text(size=20)) +
  # ggtitle("B")+
  # scale_x_discrete(name ="Distance between isolates (km)",
  # limits=c(matbind$distance_range[1:4]))+
  # scale_y_log10(labels = function(x) format(x, scientific = FALSE))+
  scale_y_continuous( trans = "log10", breaks = c(0.10,1,6),labels=c("<0.1","1.00",">6.00"),limits=c(0.1,6)) +
  # scale_y_continuous( trans = "log10", breaks = c(0.2,1,6),labels=c("<0.2","1.00",">6.00"),limits=c(0.19,6.2)) +
  
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  ylab("Risk Ratio") + 
  xlab("tMRCA (years)")+
  theme(aspect.ratio=14/15) +
  scale_color_manual(limits = c("Within Province","500-1000","Other Africa","Outside Africa"), values = c('#00539CFF','#ED2B33FF','#97BC62FF','purple')) + 
  theme(panel.grid.major =  element_line(colour="lightgrey", size=0.1), panel.grid.minor =  element_line(colour="lightgrey", size=0.1),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  theme(axis.text.y = element_text(size= 15) ) +
  xlim(0,75)
homog <- plotwithin + facet_wrap( ~ distance_range_f, nrow =1) +
  theme(strip.text.x = element_blank())



mat_sub_genomic.tmp <- subset(mat_sub_genomic, mat_sub_genomic$distance_range=="Within Province")
mat_sub_genomic.tmp$upperCI[which(mat_sub_genomic.tmp$upperCI>6)]<-6
plotwithin <- ggplot(data = mat_sub_genomic.tmp, aes( x = medMRCA, y = RR)) + 
  geom_line(color='#00539CFF')+
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI), alpha = 0.1, color = NA, fill = "#00539CFF" ) +
  theme(axis.text.x = element_text(angle = 90, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        title = element_text(size=20)) +
  scale_y_continuous( trans = "log10", breaks = c(0.10,1,6),labels=c("<0.1","1.00",">6.00"),limits=c(0.1,6)) +

  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  ylab("Risk Ratio") + 
  xlab("tMRCA (years)")+
  theme(aspect.ratio=14/15) +
  scale_color_manual(limits = c("Within Province","500-1000","Other Africa","Outside Africa"), values = c('#00539CFF','#ED2B33FF','#97BC62FF','purple')) + 
  theme(panel.grid.major =  element_line(colour="lightgrey", size=0.1), panel.grid.minor =  element_line(colour="lightgrey", size=0.1),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  theme(axis.text.y = element_text(size= 15) ) +
  xlim(0,75)
