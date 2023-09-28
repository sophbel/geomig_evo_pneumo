## Create dat.in metadata and read in all the trees  
library(ape)
library(data.table)
# ## #Read in population data for each province
load("./modelinput_data/pop_municipality.2017LS.RData")
load("./modelinput_data/cdr.mat.one.RData")
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
# ## #Read in metadata
# GPS <- read.table(file = "/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/GPS/combined_GPSPneumoDB_clean.csv", sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
# ## #Read in trees
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC2/bactdating/bd_GPSC2_1e6")
# res2 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC5/bactdating/bd_GPSC5_1e7")
# res5 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC14/bactdating/bd_GPSC14_1e7")
# res14 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC17/bactdating/bd_GPSC17_1e7")
# res17 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC13/bactdating/bd_GPSC13_1e7")
# res13 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC10/bactdating/bd_GPSC10_1e7")
# res10 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC68/bactdating/bd_GPSC68_1e8")
# res68 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC79/bactdating/sero4_14/bd_GPSC79_1e7")
# res79 <- resbd$tree
# resbd <- readRDS("/Users/sb62/Documents/Migration/BactDating/GPSC1/bactdating/bd_GPSC1_1e7")
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
# GPS_GPSC_everything <- subset(GPS, GPS$Lane_Id %in% lanes)
#  
# ##Input missing latitudes and longitudes
# # #Latitude
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "USA"] <- 37.0902
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "China"] <- 35.8617
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Egypt"] <- 26.8206 
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Thailand"] <- 15.87
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "The Gambia"] <- 13.4432
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "UK"] <- 3.4360
# 
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Austria"] <- 47.5162
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Germany"] <- 51.1657
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Iceland"] <- 64.9631
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Korea"] <- 35.9078
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Malawi"] <- 13.2543
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Malaysia"] <- 4.2105
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Netherlands"] <- 52.1326
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Portugal"] <- 39.3999
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Taiwan"] <- 23.6978
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Turkey"] <- 38.9637
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "Vietnam"] <- 14.0583
# 
# # #Longitude
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "USA"] <- -95.7129
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "China"] <- 104.195
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Egypt"] <- 30.8025
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Thailand"] <- 100.993
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "The Gambia"] <- -15.3101 
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "UK"] <- 55.3781
# 
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Austria"] <- 14.5501
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Germany"] <- 10.4515
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Iceland"] <- 19.0208
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Korea"] <- 127.7669
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Malawi"] <- 34.3015
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Malaysia"] <- 101.9758
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Netherlands"] <- 5.2913
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Portugal"] <- 8.2245
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Taiwan"] <- 120.9605
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Turkey"] <- 35.2433
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "Vietnam"] <- 108.2772
# ##Sort out the regions
# # # Subset South Africa Data and ensure region names are consistent and correct. Additionally include latitude and longitude
# GPS_GPSC_everything$Region <- as.character(GPS_GPSC_everything$Region)
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "Kwazulu Natal"] <- "KwaZulu-Natal"
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "KwaZulu Natal"] <- "KwaZulu-Natal"
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "Soweto"] <- "Gauteng"
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "Agincourt"] <- "Mpumalanga"
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "North West"] <- "North West"
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Country == "South Africa" & GPS_GPSC_everything$Region == ""] <- "South Africa"
# #Remove genomes with no region in South Africa.
# GPS_GPSC_everything <- subset(GPS_GPSC_everything, GPS_GPSC_everything$Region != "South Africa")
# # # ##LOCATION SOUTH AFRICA---------------------------------------------------------------------------
# # # #Specify that it is south africa so that it assigns the correct lat and long
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "North West"] <- "North West, South Africa"
# # # ##Province Latitude and longitude
# city <- names(table(GPS_GPSC_everything$Region[GPS_GPSC_everything$Country=="South Africa"]))
# #Geocoding
# # #input is a list of country in a csv file
# require(tmaptools)
# # #output file
# Geocoding_output <- data.frame("Region" =0, "Latitude_New"=0, "Longitude_New"=0)
# # #Geocoding
# for (i in city) {
#   longitude <- geocode_OSM(i)$coords[1]
#   latitude <- geocode_OSM(i)$coords[2]
#   Geocoding_output <- rbind(Geocoding_output, c(i,latitude, longitude))
# }
# Geocoding_output <- subset(Geocoding_output, Geocoding_output$Longitude_New != 0)
# GPS_GEO <- left_join(GPS_GPSC_everything, Geocoding_output, by = "Region")
# # # #Change it back to the name of the province once latitude is assigned
# GPS_GPSC_everything <- GPS_GEO
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Region == "North West, South Africa"] <- "North West"
# GPS_GPSC_everything$Latitude[GPS_GPSC_everything$Country == "South Africa"] <- GPS_GPSC_everything$Latitude_New[GPS_GPSC_everything$Country == "South Africa"]
# GPS_GPSC_everything$Longitude[GPS_GPSC_everything$Country == "South Africa"] <- GPS_GPSC_everything$Longitude_New[GPS_GPSC_everything$Country == "South Africa"]
# GPS_GPSC_everything$Region[GPS_GPSC_everything$Country != "South Africa"] <- GPS_GPSC_everything$Country[GPS_GPSC_everything$Country != "South Africa"]
#  
# #Drop tips
# lanes_tree <- subset(lanes, lanes  %in% GPS_GPSC_everything$Lane_Id)
# res_overall <- keep.tip(res_everything, lanes_tree)

# #Reorder tree based on order of tip labels
# order <- res_overall$tip.label 
# GPS_GPSC_everything <- GPS_GPSC_everything %>%
#   dplyr::slice(match(order, Lane_Id))
# GPS_GPSC_overall <- subset(GPS_GPSC_everything, GPS_GPSC_everything$Country == "South Africa"& GPS_GPSC_everything$Col_Year >= 2000)

# save(GPS_GPSC_overall, file = "GPS_GPSC_overall.SA.RData")
# save(GPS_GPSC_everything, file = "GPS_GPSC_everything.RData")
load("./data/gps_metadata/GPS_GPSC_overall.SA.RData")
load("./data/gps_metadata/GPS_GPSC_everything.RData")

## Create dat.in dataframe
gubgpscs <- c(1,2,5,10,13,14,17,68,79)
dat.GPSC <- output.totEvoTime.SA<-output.typesSA <- outputTipDatesSA <- TimeMRCA.SA <- output.cellindexSA <- output.meanGens.SA <- list()
for (s in gubgpscs) { 
   ##read in all GPSCs
   resbd <- readRDS(paste0("./data/phylogenies/bd_GPSC",s))
   res17 <- resbd$tree
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

save(dat.tmp.allser,file="./modelinput_data/dat.tmp.allser.RData")
save(pop_2019, file="./modelinput_data/pop_2019.RData")
