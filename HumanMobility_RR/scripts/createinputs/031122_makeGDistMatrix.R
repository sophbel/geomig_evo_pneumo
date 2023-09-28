library(ape)
library(lubridate)
library(dplyr)
load("./data/gps_metadata/GPS_GPSC_overall.SA.RData")
gubgpscs <- c(1,2,5,10,13,14,17,68,79)
tMRCAs <- gDists <- list()
for (s in gubgpscs) { 
###read in all GPSCs
resbd <- readRDS(paste0("./data/phylogenies/bd_GPSC",s))
res17 <- resbd$tree
lanes <- res17$tip.label
lanes_tree <- subset(lanes, lanes  %in% GPS_GPSC_overall$Lane_Id)
trees <- keep.tip(res17, lanes_tree)
###Reorder tree based on order of tip labels
GPS_GPSC17 <- subset(GPS_GPSC_overall, GPS_GPSC_overall$Lane_Id %in% lanes)
order <- trees$tip.label 
GPS_GPSC17 <- GPS_GPSC17 %>%
  dplyr::slice(match(order, Lane_Id))
ntrees<-length(trees$tip.label)
labels<- trees$tip.label
dates<-GPS_GPSC17$Col_time
provs<-GPS_GPSC17$Region
dat<-cbind(dates,provs,labels)

###Create empty arrays to fill
tmrcaGPSC<-sim.matsGPSC<-matrix(NaN,nrow=length(dates),ncol=length(dates))

  dist.mat<-cophenetic.phylo(trees)
  pat.dist.matGPSC<-matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))
  
  diag(pat.dist.matGPSC)<-NA
  mrca<-(outer(dates,dates,"+")-pat.dist.matGPSC)/2
  
  tmrcaGPSC<- -sweep(mrca,1,dates,"-") # Time to mrca from each tip.
  tMRCAs[[which(gubgpscs==s)]] <- tmrcaGPSC
  
attr(sim.matsGPSC,"data")<-dat

}
names(tMRCAs) <- gubgpscs

# setwd("./modelinput_data/")
# save(tMRCAs,file="tMRCAs.RData")

