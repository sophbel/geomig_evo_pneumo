## Facebook Movement Data converted to mobility between provinces for South Africa and normalized to number of facebook users per province----
'%notin%' <- Negate('%in%')
library(stringi)
library(maptools)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

### to load from here
load("./data/facebook/mvment_SA.provinces.RData")
load("./data/landscan2017/landscan_populations.RData")
mvment_SA.provinces.ZA$n_baseline_mean <- mvment_SA.provinces.ZA$n_baseline/17

# ## create mobility matrix Province----
SA_pair_BL.raw <- xtabs(n_baseline~start_province + end_province, mvment_SA.provinces.ZA)
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)

## Normalized by population
pop.prov <- with(landscan_populations, tapply(Population, list(NAME_1), FUN=sum)) ## Sum so can do it for all provinces
SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.prov[rownames(SA_pair_BL.raw.one)]  ) ) ) ) 
tmp.pop <- rowSums(SA_pair_BL.tmp)
## Forced to one
SA_pair_BL.normal.one <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )
SA_pair_BL.normal.pop <- apply(SA_pair_BL.tmp, 2, function(x) x/tmp.pop  )
probmob.Prov <- SA_pair_BL.normal.one
probmob.Prov.RAW <- SA_pair_BL.tmp
###
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




## Save file ----
cdr.mat.town.one <- probmob.Munic
cdr.mat.one <- probmob.Prov
# save(cdr.mat, file ="cdr.mat.RData")
# save(cdr.mat.town,file ="cdr.mat.town.RData")
save(cdr.mat.one, file ="./modelinput_data/cdr.mat.one.RData")
save(cdr.mat.town.one,file ="./modelinput_data/cdr.mat.town.one.RData")
# mvmentSA_surroundingcountries=mvmentSA
# save(mvmentSA_surroundingcountries, file="mvmentSA_surroundingcountries.RData")
pop2019.town<-pop.munic
save(pop2019.town,file="./modelinput_data/pop_municipality.2017LS.RData")
save(hm.prov,file="./data/facebook/mobility.plotProvince.RData")
save(hm.munic,file="./data/facebook/mobility.plotMunic.RData")


