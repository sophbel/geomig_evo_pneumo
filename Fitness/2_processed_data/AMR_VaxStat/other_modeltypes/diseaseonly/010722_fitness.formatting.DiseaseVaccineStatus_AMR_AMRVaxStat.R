library(stringr)
library(data.table)
load("/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/GPS_SA.amr.RData")
GPS_SA.amr<-subset(GPS_SA.amr,GPS_SA.amr$Manifest_type!="Carriage")
GPS_SA.amr$Region<-str_to_upper(GPS_SA.amr$Region)
vt.type <- c("NVT","VT")
amr.type <- c(0,1)
total.type <- c("nvt_s","nvt_r","vt_s","vt_r")

nb_vaccinetypes <- length(vt.type)
nb_amrtypes <- length(amr.type)
nb_types <- length(total.type)

nb_years <- length(unique(GPS_SA.amr$Year_collection))
nb_regions <- length(unique(GPS_SA.amr$Region))
years <- seq(2000,2014,1)
# years<-seq(2009,2013,1)
# provs=c("Gauteng","Mpumalanga")
provs = c("Eastern Cape", "Free State", "Gauteng", "KwaZulu Natal", "Limpopo", "Mpumalanga", "North West", "Northern Cape", "Western Cape")
provs<-str_to_upper(provs)
##########################
#####VACCINE STATUS####
##########################
clade_number_array = array(0, dim = c(nb_vaccinetypes,nb_years,  nb_regions))
clade_list<-list()
for (i in 1:nb_regions) {
  clade_number_array[,,i] <- matrix(ncol=nb_years,nrow=nb_vaccinetypes)
  clade_list[[i]] <- matrix(ncol=nb_years,nrow=nb_vaccinetypes)
  rownames(clade_list[[i]] ) <- vt.type
  # rownames(clade_list[[i]] ) <- c("NVT","PCV7","PCV13")
  colnames(clade_list[[i]] ) <- years
  # clade_number_array[,,i] <- tmp.mat
  tmp.df <- subset(GPS_SA.amr, GPS_SA.amr$Region==provs[i])
  vt.year <- table(tmp.df$Vaccine_Status, tmp.df$Year_collection)
  # vt.year <- table(tmp.df$In_Silico_Serotype, tmp.df$Col_Year)
  for (y in years) {
    for (v in vt.type) {
      clade_number_array[which(vt.type==v),which(years==y),i]  <- if (any(rownames(vt.year) == v) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == v ),which(colnames(vt.year) == y ) ] )   else  (0) 
      clade_list[[i]][which(vt.type==v),which(years==y)]  <- if (any(rownames(vt.year) == v) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == v ),which(colnames(vt.year) == y ) ] )   else  (0) 
      
    }
  }
  
}
names(clade_list)<-provs
clade_sero_list <- clade_list
setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/VaccineStatus/")
save(clade_list,file="clade_list.RData")
save(clade_sero_list,file="clade_sero_list.RData")
save(clade_number_array, file="clade_numer_array.RData")
####################################################


##########################
#####AMR####
#########################
nb_amrtypes <- length(amr.type)
clade_number_array = array(0, dim = c(nb_amrtypes,nb_years,  nb_regions))
clade_list<-list()
for (i in 1:nb_regions) {
  clade_number_array[,,i] <- matrix(ncol=nb_years,nrow=nb_amrtypes)
  clade_list[[i]] <- matrix(ncol=nb_years,nrow=nb_amrtypes)
  rownames(clade_list[[i]] ) <- amr.type
  # rownames(clade_list[[i]] ) <- c("NVT","PCV7","PCV13")
  colnames(clade_list[[i]] ) <- years
  # clade_number_array[,,i] <- tmp.mat
  tmp.df <- subset(GPS_SA.amr, GPS_SA.amr$Region==provs[i])
  vt.year <- table(tmp.df$Penicillin_WGS, tmp.df$Year_collection)
  # vt.year <- table(tmp.df$In_Silico_Serotype, tmp.df$Col_Year)
  for (y in years) {
    for (v in amr.type) {
      clade_number_array[which(amr.type==v),which(years==y),i]  <- if (any(rownames(vt.year) == v) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == v ),which(colnames(vt.year) == y ) ] )   else  (0) 
      clade_list[[i]][which(amr.type==v),which(years==y)]  <- if (any(rownames(vt.year) == v) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == v ),which(colnames(vt.year) == y ) ] )   else  (0) 
      
    }
  }
  
}
names(clade_list)<-provs
clade_sero_list <- clade_list
setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/AMR/")
save(clade_list,file="clade_list.RData")
save(clade_sero_list,file="clade_sero_list.RData")
save(clade_number_array, file="clade_numer_array.RData")
####################################################


##########################
#####AMR_VaxStat####
#########################
clade_number_array = array(0, dim = c(nb_types,nb_years,  nb_regions))
attr(clade_number_array,"vaxstat")<-c(rep(vt.type[1],2),rep(vt.type[2],2))
attr(clade_number_array,"amr")<-c(rep(amr.type,2))
attr(clade_number_array,"totaltype")<-total.type

clade_list<-list()
for (i in 1:nb_regions) {
  # clade_number_array[,,i] <- matrix(ncol=nb_years,nrow=nb_types)
  clade_list[[i]] <- matrix(ncol=nb_years,nrow=nb_types)
  rownames(clade_list[[i]] ) <- total.type
  colnames(clade_list[[i]] ) <- years
  attr(clade_list[[i]],"vaxstat")<-c(rep(vt.type[1],2),rep(vt.type[2],2))
  attr(clade_list[[i]],"amr")<-c(rep(amr.type,2))
  attr(clade_list[[i]],"totaltype")<-total.type
  for(v in vt.type){
    tmp.df1<-subset(GPS_SA.amr,GPS_SA.amr$Vaccine_Status==v)
    for(m in amr.type){
      tmp.df <- subset(tmp.df1, tmp.df1$Region==provs[i]& tmp.df1$Penicillin_WGS==m)
      vt.year <- table(tmp.df$Penicillin_WGS, tmp.df$Year_collection)
      for(y in years){
        clade_number_array[which(attr(clade_number_array,"vaxstat")==v&attr(clade_number_array,"amr")==m),which(years==y),i]<-if (any(rownames(vt.year) == m) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == m ),which(colnames(vt.year) == y ) ] )   else  (0)
        clade_list[[i]][which(attr(clade_number_array,"vaxstat")==v&attr(clade_number_array,"amr")==m),which(years==y)]  <- if (any(rownames(vt.year) == m) & any(colnames(vt.year) ==y) )  ( vt.year[which(rownames(vt.year) == m ),which(colnames(vt.year) == y ) ] )   else  (0)
      }
    }
    
  }
}


names(clade_list)<-provs
clade_sero_list <- clade_list
setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/Fitness/Data/AMR_VaxStat/diseaseonly/Overall")
save(clade_list,file="clade_list.RData")
save(clade_sero_list,file="clade_sero_list.RData")
save(clade_number_array, file="clade_numer_array.RData")





