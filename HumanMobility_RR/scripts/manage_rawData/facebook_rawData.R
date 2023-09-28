## Facebook Movement Data converted to mobility between provinces for South Africa and normalized to number of facebook users per province----
'%notin%' <- Negate('%in%')
library(stringi)
library(maptools)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
# LOAD ALL FILES IN FROM FACEBOOK MOBILITY DATA=------
# All FB mobility data from beginning of Collection in March 2020 through August 2021 excluding April 23rd to June 19th when there doesn't appear to be any data.
#
times <- c("0000","0800", "1600")
dates <- stri_pad_left(str=seq(1,25,1), 2, pad="0")
tmp <- list()
for (f in dates) {
  for (j in times) {
    tmp[[(which(dates == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_03_2021/Coronavirus_South_Africa_Jan_25_2020_Movement_between_Administrative_Regions__2021-03-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }
}
mvmentSA.tmp2 <- rbindlist(tmp)
dfsplit <- str_split_fixed(mvmentSA.tmp2$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentSA <- cbind(mvmentSA.tmp2, dfsplit)
mvmentSA$start_polygon_name <- str_to_lower(mvmentSA$start_polygon_name)
mvmentSA$end_polygon_name <- str_to_lower(mvmentSA$end_polygon_name)
### March to April 2020 Movement Data Facebook
times <- c("0000","0800", "1600")
dates <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
tmp <- list()
for (f in dates) {
  for (j in times) {
    tmp[[(which(dates == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_03_2020/1342510319081283_2020-03-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }
}
mvmentMAR <- rbindlist(tmp)
times <- c("0000","0800", "1600")
dates <- stri_pad_left(str=seq(1,23,1), 2, pad="0")
tmp2 <- list()
for (f in dates) {
  for (j in times) {
    tmp2[[(which(dates == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_03_2020/1342510319081283_2020-04-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }
}
mvmentAPR <- rbindlist(tmp2)
spring2020 <- rbind(mvmentMAR, mvmentAPR)
spring2020 <- subset(spring2020, spring2020$level == "LEVEL3" & !is.na(spring2020$n_baseline))
dfsplit <- str_split_fixed(spring2020$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
spring2020 <- cbind(spring2020, dfsplit)
spring2020$start_polygon_name <- str_to_lower(spring2020$start_polygon_name)
spring2020$end_polygon_name <- str_to_lower(spring2020$end_polygon_name)
spring2020 <- spring2020[,-c("start_quadkey","end_quadkey")]
### June to July 2020 Movement Data Facebook
times <- c("0000","0800", "1600")
dates <- stri_pad_left(str=seq(19,30,1), 2, pad="0")
tmp <- list()
for (f in dates) {
  for (j in times) {
    tmp[[(which(dates == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_06_2020/1342510319081283_2020-06-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }
}
mvmentJune <- rbindlist(tmp)
times <- c("0000","0800", "1600")
dates <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
tmp2 <- list()
for (f in dates) {
  for (j in times) {
    tmp2[[(which(dates == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_06_2020/1342510319081283_2020-07-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }
}
mvmentJuly <- rbindlist(tmp2)
summer2020 <- rbind(mvmentJune, mvmentJuly)
summer2020 <- subset(summer2020, summer2020$level == "LEVEL3" & !is.na(summer2020$n_baseline))
dfsplit <- str_split_fixed(summer2020$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
summer2020 <- cbind(summer2020, dfsplit)
summer2020$start_polygon_name <- str_to_lower(summer2020$start_polygon_name)
summer2020$end_polygon_name <- str_to_lower(summer2020$end_polygon_name)
summer2020 <- summer2020[,-c("start_quadkey","end_quadkey")]
### August to November 2020 Movement Data Facebook
months <- c("08","09","10","11")
times <- c("0000","0800", "1600")
dates <- list()
dates[[1]] <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
dates[[2]] <- stri_pad_left(str=seq(1,30,1), 2, pad="0")
dates[[3]]<- stri_pad_left(str=seq(1,31,1), 2, pad="0")
dates[[4]] <- stri_pad_left(str=seq(1,8,1), 2, pad="0")

tmp <- list()
tmp.month <- list()
for (m in 1:length(months)) {
  for (f in dates[[m]]) {
    for (j in times) {
      tmp[[(which(dates[[m]] == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_08_2020/1342510319081283_2020-",months[m],"-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    }
  }
  tmp.month[[m]] <- rbindlist(tmp)
}

mvmentAUGNOV <- rbindlist(tmp.month)
mvmentAUGNOV <- subset(mvmentAUGNOV, mvmentAUGNOV$level == "LEVEL3" & !is.na(mvmentAUGNOV$n_baseline))
dfsplit <- str_split_fixed(mvmentAUGNOV$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentAUGNOV <- cbind(mvmentAUGNOV, dfsplit)
mvmentAUGNOV$start_polygon_name <- str_to_lower(mvmentAUGNOV$start_polygon_name)
mvmentAUGNOV$end_polygon_name <- str_to_lower(mvmentAUGNOV$end_polygon_name)
mvmentAUGNOV <- mvmentAUGNOV[,-c("start_quadkey","end_quadkey")]

### November 2020 to december 2020 Movement Data Facebook
months <- c("11","12")
times <- c("0000","0800", "1600")
dates <- list()
dates[[1]] <- stri_pad_left(str=seq(9,30,1), 2, pad="0")
dates[[2]] <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
tmp <- list()
tmp.month <- list()
for (m in 1:length(months)) {
  for (f in dates[[m]]) {
    for (j in times) {
      tmp[[(which(dates[[m]] == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_11_2020/1342510319081283_2020-",months[m],"-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    }
  }
  tmp.month[[m]] <- rbindlist(tmp)
}

mvmentNOVDEC <- rbindlist(tmp.month)
mvmentNOVDEC <- subset(mvmentNOVDEC, mvmentNOVDEC$level == "LEVEL3" & !is.na(mvmentNOVDEC$n_baseline))
dfsplit <- str_split_fixed(mvmentNOVDEC$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentNOVDEC <- cbind(mvmentNOVDEC, dfsplit)
mvmentNOVDEC$start_polygon_name <- str_to_lower(mvmentNOVDEC$start_polygon_name)
mvmentNOVDEC$end_polygon_name <- str_to_lower(mvmentNOVDEC$end_polygon_name)
mvmentNOVDEC <- mvmentNOVDEC[,-c("start_quadkey","end_quadkey")]
### Janury 2021 to Feruary 2021 Movement Data Facebook
months <- c("01","02")
times <- c("0000","0800", "1600")
dates <- list()
dates[[1]] <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
dates[[2]] <- stri_pad_left(str=seq(1,16,1), 2, pad="0")
tmp <- list()
tmp.month <- list()
for (m in 1:length(months)) {
  for (f in dates[[m]]) {
    for (j in times) {
      tmp[[(which(dates[[m]] == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_11_2020/1342510319081283_2021-",months[m],"-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    }
  }
  tmp.month[[m]] <- rbindlist(tmp)
}

mvmentJANFEB <- rbindlist(tmp.month)
mvmentJANFEB <- subset(mvmentJANFEB, mvmentJANFEB$level == "LEVEL3" & !is.na(mvmentJANFEB$n_baseline))
dfsplit <- str_split_fixed(mvmentJANFEB$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentJANFEB <- cbind(mvmentJANFEB, dfsplit)
mvmentJANFEB$start_polygon_name <- str_to_lower(mvmentJANFEB$start_polygon_name)
mvmentJANFEB$end_polygon_name <- str_to_lower(mvmentJANFEB$end_polygon_name)
mvmentJANFEB <- mvmentJANFEB[,-c("start_quadkey","end_quadkey")]
### February 2021 to May 2021 Movement Data Facebook
months <- c("02","03","04","05")
times <- c("0000","0800", "1600")
dates <- list()
dates[[1]] <- stri_pad_left(str=seq(17,28,1), 2, pad="0")
dates[[2]] <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
dates[[3]] <- stri_pad_left(str=seq(1,30,1), 2, pad="0")
dates[[4]] <- stri_pad_left(str=seq(1,27,1), 2, pad="0")
tmp <- list()
tmp.month <- list()
for (m in 1:length(months)) {
  for (f in dates[[m]]) {
    for (j in times) {
      tmp[[(which(dates[[m]] == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_02_2021/1342510319081283_2021-",months[m],"-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    }
  }
  tmp.month[[m]] <- rbindlist(tmp)
}

mvmentFEBMAY <- rbindlist(tmp.month)
mvmentFEBMAY <- subset(mvmentFEBMAY, mvmentFEBMAY$level == "LEVEL3" & !is.na(mvmentFEBMAY$n_baseline))
dfsplit <- str_split_fixed(mvmentFEBMAY$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentFEBMAY <- cbind(mvmentFEBMAY, dfsplit)
mvmentFEBMAY$start_polygon_name <- str_to_lower(mvmentFEBMAY$start_polygon_name)
mvmentFEBMAY$end_polygon_name <- str_to_lower(mvmentFEBMAY$end_polygon_name)
mvmentFEBMAY <- mvmentFEBMAY[,-c("start_quadkey","end_quadkey")]
### May 2021 to August 2021 Movement Data Facebook
months <- c("05","06","07","08")
times <- c("0000","0800", "1600")
dates <- list()
dates[[1]] <- stri_pad_left(str=seq(28,31,1), 2, pad="0")
dates[[2]] <- stri_pad_left(str=seq(1,30,1), 2, pad="0")
dates[[3]] <- stri_pad_left(str=seq(1,31,1), 2, pad="0")
dates[[4]] <- stri_pad_left(str=seq(1,9,1), 2, pad="0")
tmp <- list()
tmp.month <- list()
for (m in 1:length(months)) {
  for (f in dates[[m]]) {
    for (j in times) {
      tmp[[(which(dates[[m]] == f ))*(which(times ==j))]] <-  read.table(paste0("/Users/sb62/Documents/Migration/SA_Migration_110422/SA_Metadata/facebook/Movement_Admin_Regions_05_2021/1342510319081283_2021-",months[m],"-",f,"_",j,".csv"),sep =  "," , header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
    }
  }
  tmp.month[[m]] <- rbindlist(tmp)
}

mvmentMAYAUG <- rbindlist(tmp.month)
mvmentMAYAUG <- subset(mvmentMAYAUG, mvmentMAYAUG$level == "LEVEL3" & !is.na(mvmentMAYAUG$n_baseline))
dfsplit <- str_split_fixed(mvmentMAYAUG$date_time, " ", 2)
colnames(dfsplit) <- c("date","time")
mvmentMAYAUG <- cbind(mvmentMAYAUG, dfsplit)
mvmentMAYAUG$start_polygon_name <- str_to_lower(mvmentMAYAUG$start_polygon_name)
mvmentMAYAUG$end_polygon_name <- str_to_lower(mvmentMAYAUG$end_polygon_name)
mvmentMAYAUG <- mvmentMAYAUG[,-c("start_quadkey","end_quadkey")]


## LOAD IN ALL THESE SAVED FILES
# load("/Users/sb62/Documents/Migration/Mobility_transmission_SA/Mobility_transmission/version2/outputs/mvmentSA_surroundingcountries.RData")
dfsplit <- str_split_fixed(mvment_SA.provinces.ZA$date, "-", 3)
# colnames(dfsplit) <- c("year","month","day")
mvment_SA.provinces.ZA <- cbind(mvment_SA.provinces.ZA, dfsplit)
###---- END LOADING FILES IN
mvmentSA <- rbind(summer2020,mvmentSA,spring2020,mvmentAUGNOV,mvmentNOVDEC,mvmentJANFEB,mvmentFEBMAY,mvmentMAYAUG)
mvmentSA.ZA <- subset(mvmentSA,mvmentSA$country == "ZA" )

# save(mvment_SA.provinces.ZA,file ="./data/facebook/mvment_SA.provinces.RData")


# ## Landscan Data ---- 
# ##Fix names to match
load("./data/landscan2017/LandScan_PopulationN.RData")
colnames(landscan_populations) <- c("GID_0","NAME_0","GID_1","NAME_1","GID_2","NAME_2","GID_3","NAME_3","Population")
landscan_populations$NAME_3 <- tolower(landscan_populations$NAME_3)

mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "breede river"] <- "langeberg"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "breede river"] <- "langeberg"
#
# mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "breede river"] <- "breede valley"
# mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "breede river"] <- "breede valley"

mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "kagisano"] <- "kagisano molopo"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "kagisano"] <- "kagisano molopo"
mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "molopo"] <- "kagisano molopo"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "molopo"] <- "kagisano molopo"
mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "nokeng tsa taemane"] <- "city of tshwane"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "nokeng tsa taemane"] <- "city of tshwane"
mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "kungwini"] <- "city of tshwane"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "kungwini"] <- "city of tshwane"
mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "kungwini"] <- "city of tshwane"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "kungwini"] <- "city of tshwane"
mvmentSA.ZA$start_polygon_name[mvmentSA.ZA$start_polygon_name == "west rand"] <- "randfontein"
mvmentSA.ZA$end_polygon_name[mvmentSA.ZA$end_polygon_name == "west rand"] <- "randfontein"

landscan_populations$NAME_3[landscan_populations$NAME_3 == "!kheis"] <- "kheis"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "//khara hais"] <- "khara hais"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "emnambithi/ladysmith"] <- "emnambithi"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "kai !garib"] <- "kai garib"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "breede river"] <- "breede valley"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "port st johns"] <- "port saint johns"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "tlokwe city council"] <- "tlokwe"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "victor khanye"] <- "delmas"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "ephraim mogale"] <- "greater marble hall"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "joe morolong"] <- "moshaweng"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "molopo"] <- "kagisano molopo"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "kagisano"] <- "kagisano molopo"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "kagisano/molopo"] <- "kagisano molopo"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "mfolozi"] <- "mbonambi"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "pixley ka seme"] <- "seme"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "tshwane"] <- "city of tshwane"
landscan_populations$NAME_3[landscan_populations$NAME_3 == "nokeng tsa taemane"] <- "city of tshwane"
if (length(subset(names(table(mvmentSA.ZA$start_polygon_name)) ,names(table(mvmentSA.ZA$start_polygon_name)) %notin%  names(table(landscan_populations$NAME_3))) ) ==0) {print("all names match")}
## Remove Langeberg
# landscan_populations <- subset(landscan_populations,landscan_populations$NAME_3 != "langeberg")
## Merge Landscan with mobility data----
## Change names to merge start location
landscan_populations_toMerge <- landscan_populations
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="NAME_3"] <- "start_polygon_name"
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="NAME_1"] <- "start_province"
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="Population"] <- "start_population"
lspop.merge <- landscan_populations_toMerge[,c("start_polygon_name","start_province","start_population")]
mvmentSA.ZA.tmp <- left_join(mvmentSA.ZA,lspop.merge,by="start_polygon_name") ## Merge by start municipality
## Change names to merge end location
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="start_polygon_name"] <- "end_polygon_name"
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="start_province"] <- "end_province"
colnames(landscan_populations_toMerge)[colnames(landscan_populations_toMerge)=="start_population"] <- "end_population"
lspop.merge <- landscan_populations_toMerge[,c("end_polygon_name","end_province","end_population")]
mvment_SA.provinces.ZA <- left_join(mvmentSA.ZA.tmp,lspop.merge,by="end_polygon_name") ## Merge by start municipality

## Make sure duplicate names are correct----
## Make sure emalahleni are correct
mvment_SA.provinces.ZA$start_province[mvment_SA.provinces.ZA$start_polygon_name=="emalahleni" & mvment_SA.provinces.ZA$start_polygon_id==870397] <- "Eastern Cape"
mvment_SA.provinces.ZA$end_province[mvment_SA.provinces.ZA$end_polygon_name=="emalahleni" & mvment_SA.provinces.ZA$end_polygon_id==870397] <- "Eastern Cape"
mvment_SA.provinces.ZA$start_province[mvment_SA.provinces.ZA$start_polygon_name=="emalahleni" & mvment_SA.provinces.ZA$start_polygon_id==870475] <- "Mpumalanga"
mvment_SA.provinces.ZA$end_province[mvment_SA.provinces.ZA$end_polygon_name=="emalahleni" & mvment_SA.provinces.ZA$end_polygon_id==870475] <- "Mpumalanga"
## Make sure Naledi are correct
mvment_SA.provinces.ZA$start_province[mvment_SA.provinces.ZA$start_polygon_name=="naledi" & mvment_SA.provinces.ZA$start_polygon_id==870551] <- "Free State"
mvment_SA.provinces.ZA$end_province[mvment_SA.provinces.ZA$end_polygon_name=="naledi" & mvment_SA.provinces.ZA$end_polygon_id==870551] <- "Free State"
mvment_SA.provinces.ZA$start_province[mvment_SA.provinces.ZA$start_polygon_name=="naledi" & mvment_SA.provinces.ZA$start_polygon_id==870526] <- "North West"
mvment_SA.provinces.ZA$end_province[mvment_SA.provinces.ZA$end_polygon_name=="naledi" & mvment_SA.provinces.ZA$end_polygon_id==870526] <- "North West"
## Start Location
mvment_SA.provinces.ZA$start_location <- paste0(mvment_SA.provinces.ZA$start_polygon_name,"_",mvment_SA.provinces.ZA$start_province)
mvment_SA.provinces.ZA$end_location <- paste0(mvment_SA.provinces.ZA$end_polygon_name,"_",mvment_SA.provinces.ZA$end_province)
# remove extraneous files
# rm(mvmentJuly,mvmentJune,mvmentSA.ZA,mvmentSA.tmp2,mvmentSA,mvmentSA.ZA.tmp)
# Create mean of baseline----
### save here
# save(mvment_SA.provinces.ZA,file ="./data/facebook/mvment_SA.provinces.RData")
# save(landscan_populations,file="./data/landscan2017/landscan_populations.RData")



