'%notin%'<-Negate('%in%')
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR")
#----- Script to extract Landscan population data -----#
library(maptools)
library(sf)
library(rgdal)
library(raster)
library(sp)
library(ggplot2)
library(rgeos)
library(data.table)
library(tidyr)
library(stringr)
################SET UP RR #############
#--- read in shapefile
Sa_shp <- (file = "./data/shapefiles/gadm36_ZAF_3.shp")
shp <- sf::st_read(Sa_shp)


load("./modelinput_data/pop_municipality.2017LS.RData") 
load("./ModelProjections/data/TranMatArray.234x234_MultTree_adj.RData")
TranMatArray.1<-TranMatArray.234x234
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
data.table(tn1)
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
###################population density calculation

matpop<-data.frame("NAME_3"=names(pop2019.town),"population"=pop2019.town)
shp <- merge(shp,matpop,by="NAME_3")

shp$area<-units::set_units(st_area(shp), km^2)
# shp$area<-st_area(shp)/(1000^2)
shp$pop_density<- shp$population/shp$area

shp.trans <- shp %>%
  st_geometry() %>%
  st_transform(crs = '+proj=aeqd ') 
shp.trans_simpl <- st_simplify(shp.trans, preserveTopology = FALSE, dTolerance = 1000)
# centered at input data for low distortion
shp.tmp <- shp
shp.tmp$geometry <- shp.trans_simpl


####sample population density
densities.tmp<-as.numeric(shp.tmp$pop_density)
names(densities.tmp)<-shp.tmp$NAME_3
densities<-densities.tmp[names(pop2019.town)]
# save(densities,file="./ModelProjections/data/densities.RData")
# load(file="./ModelProjections/data/densities.RData")

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
###############################################################################################
##################### Risk ratio of being at each municipality at 1 year#########
   
  nboot=1000000
   ngens=10
   noPop=FALSE
   df.whereat<-matrix(nrow=1,ncol=nboot)
   for( boot in 2:nboot) {
     if(noPop==TRUE){specificStart <- sample(234,1)}else
     {specificStart <- sample(234,1,prob=c(pop2019.town)) }
     df.whereat[1,1] <- wherenext<- specificStart
     # for (gen in 2:ngens){
       df.whereat[1,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,10])
       # df.whereat[1,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,104])
     # }
     print(boot)
   }
   
   vecloc <- vector(mode="numeric",length=234)
   df.whereat.tab<-table(df.whereat)
   for (i in 1:234) { 
      vecloc[i] <- (df.whereat.tab[i]/nboot)/mean(df.whereat.tab/nboot)
     # vecloc[i] <- (table(df.whereat)[i]/nboot)/mean(table(df.whereat)/nboot)
     print(i)
   }
   
   vecloc <- data.table(vecloc)
   vecloc$NAME_3 <- tn1
   colnames(vecloc) <- c("Risk_1year","NAME_3")
   vecloc.pop <- cbind(vecloc,pop2019.town)
   
  
   
   shp_simple <- merge(shp.tmp,vecloc,by="NAME_3")
   
# save(j_lalo,file="./ModelProjections/data/j_lalo.RData")
# save(tshw_lalo,file="./ModelProjections/data/tshw_lalo.RData")
# save(ekur_lalo,file="./ModelProjections/data/ekur_lalo.RData")
# save(ct_lalo,file="./ModelProjections/data/ct_lalo.RData")
# save(kzn_lalo,file="./ModelProjections/data/kzn_lalo.RData")
# save(nmb_lalo,file="./ModelProjections/data/nmb_lalo.RData")

# if(noPop==TRUE){shp_simple_noPop<-shp_simple
# save(shp_simple_noPop,file="./ModelProjections/data/shp_simple_noPop.RData")}else{save(shp_simple,file="./ModelProjections/data/shp_simple.RData")
# }

   ##########plot normal
risk1Year <- ggplot(data=shp_simple)+
     geom_sf(data= shp_simple,aes(fill=log(Risk_1year)), lwd = 0) +
     theme(legend.position = "none")+
  
  ### Population >1 million
  geom_sf(data=j_lalo,size=6,shape=18,color="black")+
  geom_sf(data=tshw_lalo,size=6,shape=18,color="black")+
  geom_sf(data=ekur_lalo,size=6,shape=18,color="black")+
  geom_sf(data=ct_lalo,size=6,shape=18,color="black")+
  geom_sf(data=kzn_lalo,size=6,shape=18,color="black")+
  geom_sf(data=nmb_lalo,size=6,shape=18,color="black")+
  # scale_shape_manual(values = shapes,breaks=c("Capitals","Population >1million"),limits=c("Capitals","Population >1million"))+
  theme_classic() +
  scale_fill_distiller( palette ="RdBu", direction = -1,breaks=c(-5,0,2.99),
                        labels=c(0.01,1,20),limits=c(min( log(shp_simple$Risk_1year)),log(max(shp_simple$Risk_1year))) )+
     theme(axis.text=element_blank(),
           plot.subtitle = element_text(color = "blue"),
           plot.caption = element_text(color = "Gray60"),
           legend.text=element_text(size=18),
           # legend.position = c(.9,.15),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
  legend.position = c(.12,.86),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +

     guides(fill = guide_colorbar(title = "Relative Risk",title.position = "top",
                                  title.theme = element_text(size = 18,
                                                             colour = "gray70",angle = 0)))

risk1Year

# ggsave("./ModelProjections/plots/RR1Year.map.Pop_adj.pdf",width = 10,height = 10)
upoverall<-vecloc.pop[which(vecloc.pop$Risk_1year>1)]
dim(upoverall)[1]
tmp<-mat.numIncRisk[1:18,]
quantile(upoverall$Risk_1year,probs=c(0.025,0.5,0.975))
mean(upoverall$Risk_1year)
mean(vecloc.pop[which(vecloc.pop$pop2019.town>3000000)]$Risk_1year)

 ##################################################################################


 
#######################################################
################RURAL VS URBAN###########
############################################




# ggsave("./ModelProjections/plots/ruralurban.risk1year.pdf",width=25,height=7)



###general populations
p<-ggplot(data=shp.tmp)+
   geom_sf(data= shp.tmp,aes(fill=as.numeric(pop_density)), lwd = 0)+
   scale_fill_viridis_c(option = "plasma", trans = "sqrt",
                        breaks=c(50,500,2000),limits=c(min(shp.tmp$pop_density),max(shp.tmp$pop_density))) +
   theme_light() +
   theme(axis.text=element_blank(),
         plot.subtitle = element_text(color = "blue"),
         plot.caption = element_text(color = "Gray60"),
         legend.text=element_text(size=20))  +
  annotation_scale()+
   labs(fill="Population Density\n(person/km^2)")
# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/RR1Year_popDensities.RData")
p


############ ############ ############ ############################################
####### Sequential simulations not RR for projections ###########
############ ############ ############ ####################################
#Calculating the distance travelled from Sequential simulations Reff=1
load("./modelinput_data/pairwise_geodist.town.RData")
diag(pairwise_geodist.town)<-NA
load(file="./ModelProjections/data/densities.RData")

for (i in 1){
  # nboot=600000
  ####overall
  nboot=1000
  ngens=120
  noPop=FALSE
  df.whereat<-matrix(nrow=ngens,ncol=nboot)
  munic_nums<-munic_dist<-matrix(nrow=ngens, ncol=nboot)
  
  for( boot in 2:nboot) {
    if(noPop==TRUE){specificStart <- sample(234,1)}else
    {specificStart <- sample(234,1,prob=c(pop2019.town)) }
    df.whereat[1,boot] <- wherenext<- specificStart
    for (gen in 2:ngens){
      df.whereat[gen,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,gen])
      munic_nums[gen,boot]<-length(table(df.whereat[1:gen,boot]))
      munic_dist[gen,boot]<-pairwise_geodist.town[specificStart,wherenext]
    }
    print(boot)
  }
  munic_df<-melt(munic_nums,id=1:3)
  colnames(munic_df)<-c("generations","boot","nmunic")
  munic_df<-data.frame(munic_df)
  #nmunic
  overall_avs<- data.table(apply(munic_nums, 1, function(x) mean(x,na.rm=T)))
  overall_avs$generations<-seq(1,ngens,1)
  overall_avs$years<-overall_avs$generations/10
  overall_avs$lower_CI<-apply(munic_nums, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs$upper_CI<-apply(munic_nums, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs$med<-apply(munic_nums, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  # saveRDS(overall_avs,file="Figure3/overall_avs.RData")
  
  # overall_avs$lower_CI<-apply(munic_nums, 1, function(x) min(x,na.rm = T))
  # overall_avs$upper_CI<-apply(munic_nums, 1, function(x)  max(x,na.rm = T))
  munic_df$years<-munic_df$generations/10
  
  ##distance
  munic_df2<-melt(munic_dist,id=1:3)
  colnames(munic_df2)<-c("generations","boot","dist")
  munic_df2<-data.frame(munic_df2)
  overall_avs2<- data.table(apply(munic_dist, 1, function(x) mean(x,na.rm=T)))
  overall_avs2$generations<-seq(1,ngens,1)
  overall_avs2$years<-overall_avs$generations/10
  overall_avs2$lower_CI<-apply(munic_dist, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs2$upper_CI<-apply(munic_dist, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs2$med<-apply(munic_dist, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  
  overall_avs2$min<-apply(munic_dist, 1, function(x) min(x,na.rm=T))
  overall_avs2$max<-apply(munic_dist, 1, function(x) max(x,na.rm=T))
  munic_df2$years<-munic_df2$generations/10
  # saveRDS(munic_df,file="./munic_df.RData")
  
  ggplot(munic_df)+
    geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
    geom_line(data=overall_avs,aes(x=years,y=V1))+
    theme_classic()+
    xlab("Years")+
    xlim(0,10)+
    ylab("Number of Municipalities")
  
  
  #### URBAN##########
  df.whereat<-matrix(nrow=ngens,ncol=nboot)
  munic_nums_urban<-munic_dist_urban<-matrix(nrow=ngens, ncol=nboot)
  samp<-which(densities>500)
  for( boot in 2:nboot) {
    if(noPop==TRUE){specificStart <- sample(samp,1)}else
    {specificStart <- sample(samp,1,prob=c(pop2019.town[samp])) }
    df.whereat[1,boot] <- wherenext<- specificStart
    for (gen in 2:ngens){
      df.whereat[gen,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,gen])
      munic_nums_urban[gen,boot]<-length(table(df.whereat[1:gen,boot]))
      munic_dist_urban[gen,boot]<-pairwise_geodist.town[specificStart,wherenext]
    }
    print(boot)
  }
  
  ##nmunic
  munic_df_urban<-melt(munic_nums_urban,id=1:3)
  colnames(munic_df_urban)<-c("generations","boot","nmunic")
  munic_df_urban<-data.frame(munic_df_urban)
  overall_avs_urban<- data.table(apply(munic_nums_urban, 1, function(x) mean(x,na.rm=T)))
  overall_avs_urban$generations<-seq(1,ngens,1)
  overall_avs_urban$years<-overall_avs$generations/10
  overall_avs_urban$lower_CI<-apply(munic_nums_urban, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs_urban$upper_CI<-apply(munic_nums_urban, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs_urban$med<-apply(munic_nums_urban, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  
  # overall_avs_urban$lower_CI<-apply(munic_nums_urban, 1, function(x) min(x,na.rm=T))
  # overall_avs_urban$upper_CI<-apply(munic_nums_urban, 1, function(x) max(x,na.rm=T))
  munic_df_urban$years<-munic_df_urban$generations/10
  
  ##distance
  munic_df_urban2<-melt(munic_dist_urban,id=1:3)
  colnames(munic_df_urban2)<-c("generations","boot","dist")
  munic_df_urban2<-data.frame(munic_df_urban2)
  overall_avs_urban2<- data.table(apply(munic_dist_urban, 1, function(x) mean(x,na.rm=T)))
  overall_avs_urban2$generations<-seq(1,ngens,1)
  overall_avs_urban2$years<-overall_avs$generations/10
  overall_avs_urban2$lower_CI<-apply(munic_dist_urban, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs_urban2$upper_CI<-apply(munic_dist_urban, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs_urban2$med<-apply(munic_dist_urban, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  
  overall_avs_urban2$min<-apply(munic_dist_urban, 1, function(x) min(x,na.rm=T))
  overall_avs_urban2$max<-apply(munic_dist_urban, 1, function(x) max(x,na.rm=T))
  munic_df_urban2$years<-munic_df_urban2$generations/10
  
  ggplot(munic_df_urban)+
    geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
    geom_line(data=overall_avs_urban,aes(x=years,y=V1))+
    theme_classic()+
    xlab("Years")+
    ylab("Number of Municipalities")
  
  
  
  
  ##################### RURAL############################
  df.whereat<-matrix(nrow=ngens,ncol=nboot)
  munic_nums_rural<-munic_dist_rural<-matrix(nrow=ngens, ncol=nboot)
  samp<-which(densities<=500)
  for( boot in 2:nboot) {
    if(noPop==TRUE){specificStart <- sample(samp,1)}else
    {specificStart <- sample(samp,1,prob=c(pop2019.town[samp])) }
    
    df.whereat[1,boot] <- wherenext<- specificStart
    for (gen in 2:ngens){
      df.whereat[gen,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,gen])
      munic_nums_rural[gen,boot]<-length(table(df.whereat[1:gen,boot]))
      munic_dist_rural[gen,boot]<-pairwise_geodist.town[specificStart,wherenext]
    }
    print(boot)
  }
  
  ### number of munic.
  munic_df_rural<-melt(munic_nums_rural,id=1:3)
  colnames(munic_df_rural)<-c("generations","boot","nmunic")
  munic_df_rural<-data.frame(munic_df_rural)
  overall_avs_rural<- data.table(apply(munic_nums_rural, 1, function(x) mean(x,na.rm=T)))
  overall_avs_rural$generations<-seq(1,ngens,1)
  overall_avs_rural$years<-overall_avs$generations/10
  overall_avs_rural$lower_CI<-apply(munic_nums_rural, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs_rural$upper_CI<-apply(munic_nums_rural, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs_rural$med<-apply(munic_nums_rural, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  
  # overall_avs_rural$lower_CI<-apply(munic_nums_rural, 1, function(x) min(x,na.rm=T))
  # overall_avs_rural$upper_CI<-apply(munic_nums_rural, 1, function(x) max(x,na.rm=T))
  munic_df_rural$years<-munic_df_rural$generations/10
  
  ### distance
  munic_df_rural2<-melt(munic_dist_rural,id=1:3)
  colnames(munic_df_rural2)<-c("generations","boot","dist")
  munic_df_rural2<-data.frame(munic_df_rural2)
  overall_avs_rural2<- data.table(apply(munic_dist_rural, 1, function(x) mean(x,na.rm=T)))
  overall_avs_rural2$generations<-seq(1,ngens,1)
  overall_avs_rural2$years<-overall_avs$generations/10
  overall_avs_rural2$lower_CI<-apply(munic_dist_rural, 1, function(x) quantile(x,probs=c(0.025),na.rm=T))
  overall_avs_rural2$upper_CI<-apply(munic_dist_rural, 1, function(x) quantile(x,probs=c(0.975),na.rm=T))
  overall_avs_rural2$med<-apply(munic_dist_rural, 1, function(x) quantile(x,probs=c(0.5),na.rm=T))
  
  
  overall_avs_rural2$min<-apply(munic_dist_rural, 1, function(x) min(x,na.rm=T) )
  overall_avs_rural2$max<-apply(munic_dist_rural, 1, function(x) max(x,na.rm=T))
  munic_df_rural2$years<-munic_df_rural2$generations/10
  
  ggplot(munic_df_rural)+
    geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
    geom_line(data=overall_avs_rural,aes(x=years,y=V1))+
    theme_classic()+
    xlab("Years")+
    ylab("Number of Municipalities")
  
  ggplot(munic_df_rural2)+
    geom_line(aes(x=years,y=dist,group=boot),alpha=0.09,color="grey")+
    geom_line(data=overall_avs_rural2,aes(x=years,y=V1))+
    theme_classic()+
    xlab("Years")+
    ylab("Number of Municipalities") 
}


#### Interrogate values for each at 1 year and 10 years for the number of municipalities and distance
rbind(
  overall_avs[overall_avs$years==1,],
  overall_avs_urban[overall_avs_urban$years==1,],
  overall_avs_rural[overall_avs_rural$years==1,])

rbind(
  overall_avs[overall_avs$generations==2,],
  overall_avs_urban[overall_avs_urban$generations==2,],
  overall_avs_rural[overall_avs_rural$generations==2,])

rbind(
  overall_avs[overall_avs$years==10,],
  overall_avs_urban[overall_avs_urban$years==10,],
  overall_avs_rural[overall_avs_rural$years==10,])


rbind(
  overall_avs2[overall_avs2$years==1,],
  overall_avs_urban2[overall_avs_urban2$years==1,],
  overall_avs_rural2[overall_avs_rural2$years==1,])

rbind(
  overall_avs2[overall_avs2$generations==2,],
  overall_avs_urban2[overall_avs_urban2$generations==2,],
  overall_avs_rural2[overall_avs_rural2$generations==2,])

rbind(
  overall_avs2[overall_avs2$years==10,],
  overall_avs_urban2[overall_avs_urban2$years==10,],
  overall_avs_rural2[overall_avs_rural2$years==10,])


ggplot()+
  # geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
  geom_line(data=overall_avs_rural,aes(x=generations,y=V1),color="green")+
  geom_line(data=overall_avs_urban,aes(x=generations,y=V1),color="orange")+
  geom_line(data=overall_avs,aes(x=generations,y=V1),color="purple", linetype="dashed")+
  
  theme_classic()+
  xlab("Years")+
  ylab("Number of Municipalities")

ggplot()+
  # geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
  geom_line(data=overall_avs_rural2,aes(x=generations,y=V1),color="green")+
  # geom_ribbon(data=overall_avs_rural2,aes(x=years,ymin=lower_CI,ymax=upper_CI),color="green",alpha=0.03,color=NA)+
  
  geom_line(data=overall_avs_urban2,aes(x=generations,y=V1),color="orange")+
  # geom_ribbon(data=overall_avs_urban2,aes(x=years,ymin=lower_CI,ymax=upper_CI),color="orange",alpha=0.03,color=NA)+
  # 
  geom_line(data=overall_avs2,aes(x=generations,y=V1),color="purple", linetype="dashed")+
  # geom_ribbon(data=overall_avs2,aes(x=years,ymin=lower_CI,ymax=upper_CI),color="purple",alpha=0.03,color=NA)+
  
  
  theme_classic()+
  xlab("Years")+
  ylab("Distance (km)")


hist(pairwise_geodist.town[which(densities<500),which(densities<500)],breaks=50)
hist(pairwise_geodist.town[which(densities>500),which(densities>500)],breaks=50)


samp<-which(densities<=500)
specificStart <- sample(samp,1,prob=c(pop2019.town[samp]))
wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,gen])
wherenext
pairwise_geodist.town[specificStart,wherenext]


##################################################################
