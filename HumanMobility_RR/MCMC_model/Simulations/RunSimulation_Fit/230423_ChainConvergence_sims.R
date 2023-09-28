library(data.table)
library(ggplot2)
###### ###### ###### ###### ###### ###### 
###### chain convergence for true down sample ######
###### ###### ###### ###### ###### ###### 
# iters=20010
iters=10000
thin=1

load(paste0("./MCMC_model/Simulations/output/true_downSamp/dat.in20.06.10000.1_adj.RData"))

load(paste0("./MCMC_model/Simulations/output/true_downSamp/dat.in.all0.06.10000.1_adj.RData"))

load(paste0("./MCMC_model/Simulations/output/true_downSamp/sim.ans.0.06.",iters,".1_adj.RData"))
# load(paste0("./MCMC_model/Simulations/output/true_downSamp_rescale/sim.ans.0.05.",iters,".1_adj_munic2.RData"))
chain1<-sim.ans
load(paste0("./MCMC_model/Simulations/output/true_downSamp/sim.ans.0.06.",iters,".2_adj.RData"))
# load(paste0("./MCMC_model/Simulations/output/true_downSamp_rescale/sim.ans.0.05.",iters,".2_adj_munic2.RData"))
chain2<-sim.ans
load(paste0("./MCMC_model/Simulations/output/true_downSamp/sim.ans.0.06.",iters,".3_adj.RData"))
# load(paste0("./MCMC_model/Simulations/output/true_downSamp_rescale/sim.ans.0.05.",iters,".3_adj_munic2.RData"))
chain3<-sim.ans




chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")

chains.grav<-rbind(chains.par1,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)



p<-ggplot(data=chains.grav)+
  geom_line(data=chains.grav,aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
  theme_bw()+
  ylab("pars")
p+facet_wrap(par~.)

###### ###### ###### ###### ###### ###### ###### 
###### chain convergence for random down sample ######
###### ###### ###### ###### ###### ###### ###### 
iters=10000
thin=1


load(paste0("./MCMC_model/Simulations/output/rand_downSamp/sim.ans.0.06.",iters,".1_adj2.RData"))
chain1<-sim.ans
load(paste0("./MCMC_model/Simulations/output/rand_downSamp/sim.ans.0.06.",iters,".2_adj2.RData"))
chain2<-sim.ans
load(paste0("./MCMC_model/Simulations/output/rand_downSamp/sim.ans.0.06.",iters,".3_adj2.RData"))
chain3<-sim.ans


chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")

chains.grav<-rbind(chains.par1,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)



p<-ggplot(data=chains.grav)+
  geom_line(data=chains.grav,aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
  theme_bw()+
  ylab("pars")
p+facet_wrap(par~.)

###### ###### ###### ###### ###### ###### 
#####wrong generation time. 2X######
###### ###### ###### ###### ###### ###### 
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
load("./MCMC_model/Simulations/output/wrongGenTime/sim.ans.0.06.10000.1_adj.RData")
chain1<-sim.ans
load("./MCMC_model/Simulations/output/wrongGenTime/sim.ans.0.06.10000.2_adj.RData")
chain2<-sim.ans

load("./MCMC_model/Simulations/output/wrongGenTime/sim.ans.0.06.10000.3_adj.RData")
chain3<-sim.ans
iters=10000
thin=1


chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")

chains.Munic9<-rbind(chains.par1,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)

p<-ggplot(data=chains.Munic9)+
  geom_line(data=chains.Munic9,aes(x=iters,y=chain1,group=par),color="red",alpha=0.6)+
  geom_line(data=chains.Munic9,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.6)+
  geom_line(data=chains.Munic9,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.6)+
  theme_bw()+
  ylab("pars")
p<-ggplot(data=subset(chains.Munic9,chains.Munic9$par=='par1'))+
  geom_line(data=subset(chains.Munic9,chains.Munic9$par=='par1'),aes(x=iters,y=chain1,group=par),color="red",alpha=0.6)+
  geom_line(data=subset(chains.Munic9,chains.Munic9$par=='par1'),aes(x=iters,y=chain2,group=par),color="blue",alpha=0.6)+
  geom_line(data=subset(chains.Munic9,chains.Munic9$par=='par1'),aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.6)+
  theme_bw()+
  ylab("pars")

p+facet_wrap(par~.) 

chains.burninRM <- subset(chains.Munic9,chains.Munic9$iters>1000)

p<-ggplot(data=subset(chains.burninRM,chains.burninRM$par=="par1"))+
  geom_line(aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
  geom_line(aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
  geom_line(aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
  theme_bw()+
  ylab("pars")

p+facet_wrap(par~.) 


###### ###### ###### ###### ###### ###### ###### 
###### chain convergence for only carriage####
###### ###### ###### ###### ###### ###### ###### 
iters=10000
thin=1


load(paste0("./MCMC_model/Simulations/output/carriageLocs/sim.ans.0.06.",iters,".1_adj2.RData"))
chain1<-sim.ans
load(paste0("./MCMC_model/Simulations/output/carriageLocs/sim.ans.0.06.",iters,".2_adj2.RData"))
chain2<-sim.ans
load(paste0("./MCMC_model/Simulations/output/carriageLocs/sim.ans.0.06.",iters,".3_adj2.RData"))
chain3<-sim.ans
# chain3<-chain2


chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")

chains.grav<-rbind(chains.par1,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)

p<-ggplot(data=chains.grav)+
   geom_line(data=chains.grav,aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
   geom_line(data=chains.grav,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
   geom_line(data=chains.grav,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
   theme_bw()+
   ylab("pars")
p+facet_wrap(par~.)




############################################
##### Changing the Generation Time threshold
###### ###### ###### ###### ###### ###### ###### 
###### chain convergence with changed gentime threshold####
###### ###### ###### ###### ###### ###### ###### 
iters=20000
thin=1
gens=c(15,30,60,90,120)
plotlist<-list()
for (f in 1:length(gens)){
load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.06.",iters,".1_adj_genthres",gens[f],".RData"))
chain1<-sim.ans
load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.06.",iters,".2_adj_genthres",gens[f],".RData"))
chain2<-sim.ans
load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.06.",iters,".3_adj_genthres",gens[f],".RData"))
chain3<-sim.ans
# chain3<-chain2
1-c(rejectionRate(chain1)[1],rejectionRate(chain2)[1],rejectionRate(chain3)[1])

chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")

chains.grav<-rbind(chains.par1,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)

p<-ggplot(data=chains.grav)+
  geom_line(data=chains.grav,aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
  geom_line(data=chains.grav,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
  ggtitle(gens[f])+
  theme_bw()+
  ylab("pars")
  plotlist[[f]]<-
  p+facet_wrap(par~.)

}
library(patchwork)
(plotlist[[1]]+plotlist[[2]]+plotlist[[3]] )/ (plotlist[[4]]+plotlist[[5]])

###### ###### ###### ###### ###### ###### ###### 
###### chain convergence with changed gentime threshold only estimating a single parameter####
###### ###### ###### ###### ###### ###### ###### 
iters=20000
thin=1
gens=c(15,30,60,90,120)
plotlist<-list()
for (f in 1:length(gens)){
  load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.07.",iters,".1_adj_genthres",gens[f],"_singPar.RData"))
  chain1<-sim.ans
  load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.07.",iters,".2_adj_genthres",gens[f],"_singPar.RData"))
  chain2<-sim.ans
  load(paste0("./MCMC_model/Simulations/output/genTimeThresh/sim.ans.0.07.",iters,".3_adj_genthres",gens[f],"_singPar.RData"))
  chain3<-sim.ans
  # chain3<-chain2
  print(1-c(rejectionRate(chain1)[1],rejectionRate(chain2)[1]))
  # sim.ans.0.07.20000.2_adj_genthres120_singPar.RData
  chains.par1<-data.table("chain1"=chain1[,1],"chain2"=chain2[,1],"chain3"=chain3[,1],"iters"=seq(1,iters,thin),"par"="par1")
  # chains.par2<-data.table("chain1"=chain1[,2],"chain2"=chain2[,2],"chain3"=chain3[,2],"iters"=seq(1,iters,thin),"par"="par2")
  # chains.par3<-data.table("chain1"=chain1[,3],"chain2"=chain2[,3],"chain3"=chain3[,3],"iters"=seq(1,iters,thin),"par"="par3")
  # chains.par4<-data.table("chain1"=chain1[,4],"chain2"=chain2[,4],"chain3"=chain3[,4],"iters"=seq(1,iters,thin),"par"="par4")
  # chains.par5<-data.table("chain1"=chain1[,5],"chain2"=chain2[,5],"chain3"=chain3[,5],"iters"=seq(1,iters,thin),"par"="par5")
  # chains.par6<-data.table("chain1"=chain1[,6],"chain2"=chain2[,6],"chain3"=chain3[,6],"iters"=seq(1,iters,thin),"par"="par6")
  # chains.par7<-data.table("chain1"=chain1[,7],"chain2"=chain2[,7],"chain3"=chain3[,7],"iters"=seq(1,iters,thin),"par"="par7")
  # chains.par8<-data.table("chain1"=chain1[,8],"chain2"=chain2[,8],"chain3"=chain3[,8],"iters"=seq(1,iters,thin),"par"="par8")
  # chains.par9<-data.table("chain1"=chain1[,9],"chain2"=chain2[,9],"chain3"=chain3[,9],"iters"=seq(1,iters,thin),"par"="par9")
  
  chains.grav<-rbind(chains.par1)#,chains.par2,chains.par3,chains.par4,chains.par5,chains.par6,chains.par7,chains.par8,chains.par9)
  
  p<-ggplot(data=chains.grav)+
    geom_line(data=chains.grav,aes(x=iters,y=chain1,group=par),color="red",alpha=0.8)+
    geom_line(data=chains.grav,aes(x=iters,y=chain2,group=par),color="blue",alpha=0.8)+
    geom_line(data=chains.grav,aes(x=iters,y=chain3,group=par),color="darkgreen",alpha=0.8)+
    ggtitle(gens[f])+
  theme_bw()+
    ylab("pars")
  plotlist[[f]]<-
  p+facet_wrap(par~.)
  
}

library(patchwork)
(plotlist[[1]]+plotlist[[2]]+plotlist[[3]])/(plotlist[[4]]+plotlist[[5]])





