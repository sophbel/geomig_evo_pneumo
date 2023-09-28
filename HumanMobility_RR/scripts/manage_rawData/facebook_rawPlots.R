########
### This script plots the radius of gyration from the raw facebook data

'%notin%' <- Negate('%in%')
library(stringi)
library(maptools)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
###########
##load human mobility data
##########

load("./data/facebook/raw_data/mvment_SA.provinces.RData")
load("./data/landscan2017/landscan_populations.RData")

# ## create mobility matrix Municipality----
SA_pair_BL.raw <- xtabs(n_baseline~start_location + end_location, mvment_SA.provinces.ZA)##BASELINE
# SA_pair_BL.raw <- xtabs(n_crisis~start_polygon_name + end_polygon_name, mvment_SA.provinces.ZA)##CRISIS
# ##INDEX IDS TO POLYGON NAMES OR CONCATENATE PROVINCE NAMES ON TO POLYGON NAMES
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)
## Normalized by population
landscan_populations$location <- paste0(landscan_populations$NAME_3,"_",landscan_populations$NAME_1)
pop.munic <- landscan_populations$Population
names(pop.munic) <- landscan_populations$location
# SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.munic[rownames(SA_pair_BL.raw.one)]  ) ) ) ) 
SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.munic[rownames(SA_pair_BL.raw.one)]  ) ) ) )

tmp.pop <- rowSums(SA_pair_BL.tmp)
## Forced to one
# SA_pair_BL.normal.pop <- apply(SA_pair_BL.tmp, 2, function(x) x/tmp.pop  )
SA_pair_BL.normal.one <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )

probmob.Munic <- SA_pair_BL.normal.one
probmob.Munic.RAW <- SA_pair_BL.raw


#################
####deal with map file
################
#--- read in shapefile
Sa_shp <- (file = "./data/shapefiles/gadm36_ZAF_3.shp")
shp <- sf::st_read(Sa_shp)


load("./modelinput_data/pop_municipality.2017LS.RData") 
# load("./ModelProjections/data/TranMatArray.234x234_adj.RData")
# TranMatArray.1<-TranMatArray.234x234
load("./modelinput_data/cdr.mat.town.one.RData")# # [Mobility_ManyMonthsSA.R] ## Probability of movement between each province normalized to carriage rates for each province 
cdr.mat.town<-cdr.mat.town.one
tn <- rownames(cdr.mat.town)
pop2019.town <-pop2019.town[tn]

tn[grep("emalah",tn)]<-c("emalahleni.EC_Eastern Cape","emalahleni.MP_Mpumalanga")
tn[grep("naled",tn)]<-c("naledi.FS_Free State","naledi.NW_Nort West")
tn1<-tn
for (i in 1:length(tn)){
  tn1[i] <- strsplit(tn,"_")[[i]][1]
}
tn1<-str_to_title(tn1)
# moshaweng is botswana 
shp$NAME_3<-as.character(shp$NAME_3)
shp$NAME_3[grep("Emalah",shp$NAME_3)]<-c("Emalahleni.ec","Emalahleni.mp")
shp$NAME_3[grep("Naled",shp$NAME_3)]<-c("Naledi.fs","Naledi.nw")
shp$NAME_3<-as.factor(shp$NAME_3)
shpnames <- shp$NAME_3
tn1[which(is.na(match(tn1,shpnames)))] <- c("City of Cape Town", "City of Johannesburg", "City of Matlosana", "City of Tshwane","Delmas","Dr JS Moroka","eDumbe","Emnambithi/Ladysmith","Greater Marble Hall",
                                            "Kagisano/Molopo","Kai !Garib","//Khara Hais","!Kheis","KwaDukuza","Maluti a Phofung","Mfolozi","Moshaweng","Port St Johns","Pixley Ka Seme","Tlokwe City Council",
                                            "uMhlathuze","uMlalazi","uMngeni","uMshwathi","UMuziwabantu","UPhongolo")
tn1[which(tn1=="Greater Marble Hall")] <- "Ephraim Mogale"
tn1[which(tn1=="Delmas")] <- "Victor Khanye"
tn1[which(tn1=="Moshaweng")] <- "Joe Morolong"
names(pop2019.town) <- tn1
save(tn1,file="./modelinput_data/name_match.RData")
###################population density calculation

matpop<-data.frame("NAME_3"=names(pop2019.town),"population"=pop2019.town)
shp <- merge(shp,matpop,by="NAME_3")

shp$area<-units::set_units(st_area(shp), km^2)
# shp$area<-st_area(shp)/(1000^2)
shp$pop_density<- shp$population/shp$area


shp.trans <- shp %>%
  st_geometry() %>%
  st_transform(crs = '+proj=aeqd') 
shp.trans_simpl <- st_simplify(shp.trans, preserveTopology = FALSE, dTolerance = 1000)

save(probmob.Munic.RAW,file="./modelinput_data/probmob.Munic.RAW.RData")
#########################
#################
###calculate radius of gyration
##############
load("./modelinput_data/pairwise_geodist.town.RData")
# diag(probmob.Munic.RAW)<-NA
rucm<-vector(mode="numeric",length=234)
for ( i in 1:nrow(probmob.Munic.RAW)){
  probvec<-probmob.Munic.RAW[i,]/sum(probmob.Munic.RAW[i,],na.rm=T)
  dests<-which(probmob.Munic.RAW[i,]!=0)
  Nu<-length(dests)
  sum_icm<-vector(mode="numeric",length=Nu)
    for(q in 1:length(dests)){
        sum_icm[q]<-probvec[dests[q]]*((pairwise_geodist.town[i,dests[q]])^2)
    }
  rucm[i]<-sqrt((sum(sum_icm)))
}


rg_df<-data.frame("NAME_3"=tn1,"radius_of_gyration"=rucm)
shp <- merge(shp,rg_df,by="NAME_3")

#########
###############top population sizes
## GT
j_lalo <- data.table(t(data.table(c(-26.2041, 28.0473))))## Joburg
colnames(j_lalo) <- c("latitude","longitude")
j_lalo <- st_as_sf(j_lalo, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")

## GT ### Not Capital 25.6051° S, 28.3929° E
tshw_lalo <- data.table(t(data.table(c(-25.6051, 28.3929))))## City of Tshwane
colnames(tshw_lalo) <- c("latitude","longitude")
tshw_lalo <- st_as_sf(tshw_lalo, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

## GT ### Not Capital 26.1777° S, 28.3462° E
ekur_lalo <- data.table(t(data.table(c(-26.1777, 28.3462))))## Ekurhuleni
colnames(ekur_lalo) <- c("latitude","longitude")
ekur_lalo <- st_as_sf(ekur_lalo, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

## CT
ct_lalo <- data.table(t(data.table(c(-33.9249, 18.4241)))) ## Cape Town
colnames(ct_lalo) <- c("latitude","longitude")
ct_lalo <- st_as_sf(ct_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
## NC
sp_lalo <- data.table(t(data.table(c(-28.7553, 24.6668)))) ## Sol Plaatjie
colnames(sp_lalo) <- c("latitude","longitude")
sp_lalo <- st_as_sf(sp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##NW
# 25.8560° S, 25.6403° E
nw_lalo <- data.table(t(data.table(c(-25.8560,25.6403))))
colnames(nw_lalo) <- c("latitude","longitude")
nw_lalo <- st_as_sf(nw_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##MP
# 25.4753° S, 30.9694° E
mp_lalo <- data.table(t(data.table(c(-25.4753,30.9694))))
colnames(mp_lalo) <- c("latitude","longitude")
mp_lalo <- st_as_sf(mp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##FS
# 29.0852° S, 26.1596° E
fs_lalo <- data.table(t(data.table(c(-29.0852,26.1596))))
colnames(fs_lalo) <- c("latitude","longitude")
fs_lalo <- st_as_sf(fs_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##LP
# 23.8962° S, 29.4486° E
lp_lalo <- data.table(t(data.table(c(-23.8962,29.4486))))
colnames(lp_lalo) <- c("latitude","longitude")
lp_lalo <- st_as_sf(lp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##EC
# 32.9344° S, 27.6435° E
ec_lalo <- data.table(t(data.table(c(-32.9344,27.6435))))
colnames(ec_lalo) <- c("latitude","longitude")
ec_lalo <- st_as_sf(ec_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")


##EC ### Port Elizabeth - not the capital Nelson Mandela Bay
# 33.7452° S, 25.5681° E
nmb_lalo <- data.table(t(data.table(c(-33.7452,25.5681))))
colnames(nmb_lalo) <- c("latitude","longitude")
nmb_lalo <- st_as_sf(nmb_lalo, coords = c("longitude", "latitude"), 
                     crs = 4326, agr = "constant")

##KZN ### Ethekwini ### Durban
# 29.8587° S, 31.0218° E
kzn_lalo <- data.table(t(data.table(c(-29.8587,31.0218))))
colnames(kzn_lalo) <- c("latitude","longitude")
kzn_lalo <- st_as_sf(kzn_lalo, coords = c("longitude", "latitude"), 
                     crs = 4326, agr = "constant")



shapes=c("Capitals"=8,"Population >1million"=18)
##################

p<-ggplot(data=shp)+
  geom_sf(data= shp,aes(fill=radius_of_gyration), lwd = 0)+
  scale_fill_viridis_c(option = "plasma", trans = "sqrt",
                       breaks=c(0,200,600),limits=c(min(shp$radius_of_gyration),600)) +
  theme_light() +
  theme(axis.text=element_blank(),
        plot.subtitle = element_text(color = "blue"),
        plot.caption = element_text(color = "Gray60"),
        legend.text=element_text(size=10))  +
  annotation_scale()+
  ### Population >1 million
  # geom_sf(data=j_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_sf(data=tshw_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_sf(data=ekur_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_sf(data=ct_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_sf(data=kzn_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_sf(data=nmb_lalo,size=6,shape=18,color="grey",alpha=0.7)+
  # geom_circle()+
  labs(fill="Radius of Gyration (km)")
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/RR1Year_popDensities.RData")
p



# 
# 
# shporder<-shp$NAME_3
# # centered at input data for low distortion
# #### reorder to be in the same order as the mobility data
# shp.trans_simpl<-shp.trans_simpl[order(match(shporder,tn1)),]
# ### calculate mean distance in km
# vec_dist[1]<-mean(st_distance(st_centroid(shp.trans_simpl[1]),st_centroid(shp.trans_simpl[dests])))/1000
# 
# ####################








####################
### merge in mobility numbers
####################
# 
# ##find the probability that each municipality is a destination overall
# destination_probstmp<-colSums(probmob.Munic.RAW)
# destination_probs<-(destination_probstmp/sum(destination_probstmp))
# ##find the probability that each municipality is an origin overall
# origin_probstmp<-rowSums(probmob.Munic.RAW)
# origin_probs<-(origin_probstmp/sum(origin_probstmp))
# ##fbind them to merge into shp.tmp
# raw_mobility<-data.table(cbind(rownames(probmob.Munic.RAW),as.numeric(origin_probs),as.numeric(destination_probs)))
# raw_mobility$V2<-as.numeric(raw_mobility$V2)
# raw_mobility$V3<-as.numeric(raw_mobility$V3)
# colnames(raw_mobility)<-c("NAME_3","Origin_prob","Destination_prob")
# raw_mobility$NAME_3<-tn1
# shp_simple <- merge(shp.tmp,raw_mobility,by="NAME_3")
# 
# ##########plot origin probability
# origin_prob <- ggplot(data=shp_simple)+
#   geom_sf(data= shp_simple,aes(fill=Destination_prob), lwd = 0) +
#   theme(legend.position = "none")+
#   theme_classic() +
#   scale_fill_distiller( palette ="RdBu", direction = -1,breaks=c(0,0.001,3),
#                         labels=c(0,0.001,3),limits=c(min( shp_simple$Destination_prob),max(shp_simple$Destination_prob)) )+
#   theme(axis.text=element_blank(),
#         plot.subtitle = element_text(color = "blue"),
#         plot.caption = element_text(color = "Gray60"),
#         legend.text=element_text(size=18),
#         # legend.position = c(.9,.15),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
#         legend.position = c(.12,.86),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
#   
#   guides(fill = guide_colorbar(title = "Origin Probability",title.position = "top",
#                                title.theme = element_text(size = 18,
#                                                           colour = "gray70",angle = 0)))
# origin_prob
# 



