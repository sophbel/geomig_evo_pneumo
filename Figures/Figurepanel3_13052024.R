#############
######FIGURE 3##########
################
library(ggplot2)
####Figure 3A ####
load("./Figures/Figure3/df.propperprov.RData")
df.propperprov<-data.table(t(df.propperprov))
colnames(df.propperprov)<-c("lowerprop","prop","upperprop","popsize")
A3<-ggplot(df.propperprov,aes(x=popsize,y=prop))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lowerprop,ymax=upperprop))+
  theme_classic()+
  xlab("Population size (x10^6)")+
  ylab("Proportion of Infections")+
  ylim(0,0.4)+
  theme(axis.text = element_text(size=20),axis.title=element_text(size=20))


summary(lm(df.propperprov$prop~df.propperprov$popsize))
####Figure 3B####
load("./Figures/Figure3/dist.quants.tmp3.posteriorsAll_multTrees_adj.RData")

'%notin%'<-Negate('%in%')
genvec1<-c("Meta","Random" )
# genvec1<-c("Individual","Meta","Random" )

indivs<-subset(dist.quants.tmp3,dist.quants.tmp3$gens%in%genvec1)
trans.gens<-subset(dist.quants.tmp3,dist.quants.tmp3$gens%notin%genvec1)
trans.gens1<-subset(trans.gens,trans.gens$gens%notin%c("Individual_infec"))
B3<-ggplot() +
  geom_point(data=indivs,aes(x=distances,y=median,color=gens),size=3,position=position_dodge(width=0.5),alpha=0.8) +

  geom_point(data=trans.gens1,aes(x=distances,y=median,shape=gens),size=3,position=position_dodge(width=0.35)) +
  geom_errorbar(data=trans.gens1,aes(ymin=lowerCI,ymax=upperCI,x=distances,shape=gens),alpha=0.6,width=0.5,position=position_dodge(width=0.35))+
  # geom_errorbar(data=indivs,aes(ymin=lowerCI,ymax=upperCI,x=distances,color=gens),alpha=0.6,width=0.5,position=position_dodge(width=0.35))+
  
  xlab("Distance \n(from starting municipality)")+
  ylab("Proportion")+
  theme(panel.grid.minor  = element_blank(),)+
  theme_classic()+
  scale_shape_manual(values=c(7,19,2),limits=c("Individual_daily",1,10), labels=c("0.5","1","10")) +
  scale_color_manual(limits=c("Meta","Random"), labels=c("Meta Data","Random Mobility"),values=c("darkblue","maroon") )+
  
  ylim(0,1)+
  theme(
    # axis.text.x = element_text(angle=90,size=rel(2)), axis.text.y=element_text(size=rel(2)),
    axis.text.x = element_text(angle=90,vjust=0.6,size=20), axis.text.y=element_text(size=20),
    axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
    legend.text = element_text(size=10),legend.title = element_blank(),
    legend.position = c(.9,.7),legend.background = element_rect(color=NA, size=0.5, linetype="solid"))+
  
  # scale_shape_manual(name="Time",values=c(19,2,7),limits=c(1,10,"Random"), labels=c("1 Generation","50 Generations","Given random \nmobility")) +
  scale_x_discrete(limits=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"),labels=c("Within Municipality" ,"10-200km","200-500km","500-1000km",">1000km"))+
  guides(fill="none")
B3

####Figure 3C Map ####
shp<-sf::st_read("./Figures/Figure3/shps/gadm36_ZAF_3.shp")
load("./Figures/Figure3/shp_simple.RData")
shp_simpleTMP<-shp_simple[,c("Risk_1year","pop_density","NAME_3")]
shp_simpleTMP<-st_drop_geometry(shp_simpleTMP)
shp <- merge(shp,shp_simpleTMP,by="NAME_3")
shp_simpleplot<-shp
##########plot Risk at 1 year
risk1Year <- ggplot(data=shp_simpleplot)+
  geom_sf(data= shp_simpleplot,aes(fill=log(Risk_1year)), lwd = 0) +
  theme(legend.position = "none")+
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


###Figure 3D####
overall_avs<-readRDS(file="./Figures/Figure3/overall_avs.RData")
munic_df<-readRDS(file="./Figures/Figure3/munic_df.RData")
D3<-ggplot(munic_df)+
  geom_line(aes(x=years,y=nmunic,group=boot),alpha=0.09,color="grey")+
  geom_line(data=overall_avs,aes(x=years,y=V1))+
  theme_classic()+
  xlab("Years")+
  xlim(0,10)+
  ylab("Number of Municipalities")



###Figure 3E####
load("./Figures/Figure3/mat.nvtpcv_adj.RData")
E3<-ggplot(subset(mat.nvtpcv,mat.nvtpcv$type=="Municipalites"),aes(x=year,y=mean)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI), size=0.8) +
  theme_classic()+
  scale_y_continuous( trans = "log10" , breaks = c(0.8,1,2),labels=c("<0.8","1.00","2.00"),limits=c(0.8,2)) +
  ylab("Relative Municipalities Visited\n(Non-Vaccine Type/ Vaccine Type")+
  xlab("")+
  geom_hline(aes(yintercept=1),linetype="dashed",color="red")+
  scale_x_discrete(limits=c(1,2),labels=c("1 yr","2 yrs"))+
  theme(axis.text=element_text(size=20),axis.title = element_text(size=20),axis.text.x=element_text(size=20,angle=45,vjust=0.6))+
  labs(shape="Years of\nTransmission")


library(patchwork)
(A3+B3)/(D3+E3)
(risk1Year)


