library(tmaptools)
library(geodist)
library(PBSmapping)
library(data.table)
library(ggplot2)
library(patchwork)

## Pairwise distances on the province level
load("./modelinput_data/cdr.mat.one.RData")# # [mobility_probability.R] ## Probability of movement between each province normalized to carriage rates for each province 
cdr.mat<-cdr.mat.one
## change North West to North West South Africa so as not to confuse geoCode
provs <- colnames(cdr.mat)
provs[provs == "North West"] <- "North West, South Africa"
require(tmaptools)
# #output file
Geocoding_output <- data.frame("Region" =0, "Latitude_New"=0, "Longitude_New"=0)
# #Geocoding to assign longitude and latitudes for each province
for (i in provs) {
  longitude <- geocode_OSM(i)$coords[1]
  latitude <- geocode_OSM(i)$coords[2]
  Geocoding_output <- rbind(Geocoding_output, c(i,latitude, longitude))
}
Geocoding_output <- subset(Geocoding_output, Geocoding_output$Region!=0)
Geocoding_output$Region[Geocoding_output$Region == "North West, South Africa"] <- "North West" # Change North West back
# #CREATE PAIRWISE DISTANCE MATRIX FOR REGIONS
# # POPULATE DISTANCE MATRIX WITH CONTINOUS DISTANCES PAIRWISE -------------
lat_vec <- Geocoding_output$Latitude_New
long_vec <- Geocoding_output$Longitude_New
region_vec <- Geocoding_output$Region
x <- cbind(long_vec,lat_vec)
y <- x
pairwise_geodist <- geodist(x,y,paired = FALSE,measure = "haversine") /1000
colnames(pairwise_geodist) <- region_vec
rownames(pairwise_geodist) <- region_vec
pairwise_geodist<- round(pairwise_geodist, 2)
# save(pairwise_geodist,file = "./modelinput_data/pairwise_geodist.RData" )

load("./data/facebook/raw_data/mvment_SA.provinces.RData")

## find unique destinations
grid.coords.final <- unique(mvment_SA.provinces.ZA[,c("end_lon","end_lat")])
## unique origins
unique_municipalities <- unique(mvment_SA.provinces.ZA[,c("start_location","start_lon", "start_lat")])
colnames(unique_municipalities)<-c("name","X","Y")
attr(unique_municipalities,"projection")<-"LL"
grid.coords.final.UTM <-convUL(unique_municipalities)
dist.matrix <- as.matrix(dist(grid.coords.final.UTM[,2:3]))
rownames(dist.matrix) <- unique_municipalities$name
colnames(dist.matrix) <- unique_municipalities$name

load("./modelinput_data/cdr.mat.town.one.RData")
cdr.mat.town<-cdr.mat.town.one
mobNames <- rownames(cdr.mat.town)
pairwise_geodist.town <- dist.matrix[-c(which(duplicated(rownames(dist.matrix)))),-c(which(duplicated(rownames(dist.matrix))))]
# save(pairwise_geodist.town,file = "./modelinput_data/pairwise_geodist.town.RData" )


### plot province distance
df_prov<-matrix(nrow=nrow(pairwise_geodist)^2,ncol=2)
df_prov[,1]<-as.vector(mode="numeric",pairwise_geodist)
df_prov<-data.table(df_prov)
df_prov[which(df_prov$V1==0),]<-NA
prov<-ggplot()+
  geom_histogram(aes(df_prov$V1),colour = "black", fill = "grey", bins = 15)+
  theme_classic()+
  theme(axis.text = element_text(size=15),axis.title = element_text(size=15))+
  xlab("Distance Between Provinces (km)\n(N=9)")

### plot munic distance
df_prov2<-matrix(nrow=nrow(pairwise_geodist.town)^2,ncol=2)
df_prov2[,1]<-as.vector(mode="numeric",pairwise_geodist.town)
df_prov2<-data.table(df_prov2)
df_prov2[which(df_prov2$V1==0),]<-NA
munic<-ggplot()+
  geom_histogram(aes(df_prov2$V1),colour = "black", fill = "grey", bins = 15)+
  theme_classic()+
  theme(axis.text = element_text(size=15),axis.title = element_text(size=15))+
  xlab("Distance Between Municipalities (km)\n(N=234)")
prov+munic
