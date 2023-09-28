##Load shapefiles
# Sa_shp <- (file = "/Users/sb62/Documents/Migration/Shapefiles/ZA_shapefiles_levels/gadm36_ZAF_3.shp")
# Sa_shp.4 <- (file = "/Users/sb62/Documents/Migration/Shapefiles/ZA_shapefiles_levels/gadm36_ZAF_4.shp")
# 
# population <- data.frame(province = c("Gauteng", "KwaZulu-Natal", "Western Cape", "Eastern Cape", "Limpopo", "Mpumalanga", "North West", "Free State", "Northern Cape"),
#                          population = c(14278669, 11074784, 6510312, 6498683, 5778442, 4444212, 3856174, 2866678, 1213996))
# shapefile <- sf::st_read(Sa_shp)
##Population 2019
#----- Script to extract Landscan population data -----#
library(maptools)
library(sf)
library(rgdal)
library(raster)

#--- read in shapefile
Sa_shp <- (file ="./data/ZA_shapefiles/gadm36_ZAF_3.shp")



shp <- sf::st_read(Sa_shp)

#--- read in raster 
pop <- raster('./data/landscan2017/raw_data/w001001.adf')
# plot(pop)

#--- check coordinate ref systems are the same
# crs(pop)
# crs(shp)
# if not same use st_transform() to transform data to same CRS, shown below:
# pop <- st_transform(shp, crs(pop))


#--- crop & mask raster to area of interest 
pop_c <- raster::crop(pop, extent(shp))
pop_c <- raster::mask(pop_c, shp)
#--- Function to extract population data from raster
extractPop <- function(pop, shp){
  
  tmpj <- raster::crop(pop, extent(shp)) # crop the raster to area of interest
  tmpj <- raster::mask(tmpj, shp) # remove values outside shapefile
  pop_n <- as.vector(tmpj) # vectorize population values
  pop_n <- na.omit(pop_n) # remove NAs
  pop_n <- pop_n[pop_n>=0] # remove any remaining negative values
  tot_pop <- round(sum(pop_n)) # sum population
  return(tot_pop)
}
#--- Dataframe to store extracted data
df <- data.frame("GID_0"=NA,"NAME_0"=NA,"GID_1"=NA,"NAME_1"=NA,"GID_2"=NA,
                 "NAME_2"=NA,"GID_3"=NA,"NAME_3"=NA,Population=NA)
shp$GID_0 <- as.character(shp$GID_0)
shp$NAME_0 <- as.character(shp$NAME_0)
shp$GID_1 <- as.character(shp$GID_1)
shp$NAME_1 <- as.character(shp$NAME_1)
shp$GID_2 <- as.character(shp$GID_2)
shp$NAME_2 <- as.character(shp$NAME_2)
shp$GID_3 <- as.character(shp$GID_3)
shp$NAME_3 <- as.character(shp$NAME_3)

#--- Loop to extract population data for each region in shapefile
for(i in 1:nrow(shp)){
  
  # subset shapefile
  shp_i <- shp[i, ]
  
  # extract pop
  pop_i <- extractPop(pop, shp_i)
  df[i,1:8] <- paste(shp_i[1,c(1:4,6:7,9:10)]) # paste names of each area to df
  df$Population[i] <- pop_i
  
  # print progress
  print(paste(i, '/', nrow(shp), sep=''))
  
}
landscan_populations <- df
save(landscan_populations,file="./data/landscan2017/LandScan_PopulationN.RData")
load("./data/landscan2017/LandScan_PopulationN.RData")
