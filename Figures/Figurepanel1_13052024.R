#############
######FIGURE 1##########
################

###load libraries
library(ggplot2)
library(data.table)

####Figure 1A###
## Figure 1A is made in Figtree and Microreact
# url{https://microreact.org/project/7wqgd2gbBBEeBLLPKonbaT-belman2024southafricapneumococcus}

###Figure 1B####
shp.tmp<-sf::st_read("./Figures/Figure1/shps/gadm36_ZAF_1.shp")
provs=c("Eastern Cape","Free State","Gauteng","KwaZulu-Natal","Limpopo","Mpumalanga","North West","Northern Cape","Western Cape")
cols = c("#1A9D75","#D95E02","#7470B3","#E7288B","#66A71E","#E7AB02","#A7761C","#616161","#00008B")
B1<-ggplot(data=shp.tmp)+
  geom_sf(data=shp.tmp,aes(fill=NAME_1) )+
  theme_classic()+
  scale_fill_manual(breaks=provs,values=cols,name="Province")+
  # ggtitle("A")+
  theme(axis.text=element_blank(),legend.text = element_text(size=15),legend.title = element_text(size=15),title = element_text(size=15),legend.position = "none")+
  # scalebar(data=shp.tmp,dist = 500,dist_unit = "km", location="bottomleft",transform = FALSE, 
           # model = "WGS84", height = 0.03, st.bottom=TRUE,
           # st.dist = 0.04,st.size = 6)+
  theme(axis.title=element_blank())


###Figure 1C####
## Collection year plot####
load("./Figures/Figure1/GPS_SA.RData")
load("./Figures/Figure1/GPS_GPSC_overall.SA.RData")
GPS_SA <- subset(GPS_SA, GPS_SA$Col_Year>1999)
GPS_SA$vax.period<-ifelse(GPS_SA$Col_Year<2009, "prePCV7","postPCV")
colors=c("South Africa"="black","Dominant GPSCs"="#BF40BF")
C1<-ggplot()+
  theme_classic()+
  geom_line(data=GPS_SA,aes(fill=..count..,x=Col_Year,color="South Africa"),size=1,stat="bin",binwidth=1)+
  geom_line(data= GPS_GPSC_overall,aes(fill=..count..,x=Col_Year,color="Dominant GPSCs"),size=1,stat="bin",binwidth=1)+
  labs(x="Collection Year",
       y="count",
       color="")+
  scale_color_manual(values=colors,breaks=c("South Africa","Dominant GPSCs"),limits=c("South Africa","Dominant GPSCs"))+
  theme(axis.text = element_text(size=15), axis.title = element_text(size=15),
        legend.text = element_text(size=15),legend.title = element_text(size=15),legend.position = "bottom",
        # plot.margin=unit(c(1,2,1,1),"cm"),
        title = element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=10)))#+ggtitle("B")

# ggsave("/Users/sophiebelman/Documents/Migration/SA_Migration_110422/Submission/Figures/Figure_panel1/collectionYear.pdf",width=5.5,height=5)


####Figure 1D####
load("./Figures/Figure1/overall_data_multiTree.RData")
load("./Figures/Figure1/true.dat.quants.Munic_multTree.RData")
load("./Figures/Figure1/dat.quants.Munic_multTree.RData")

##########Visualize Fit####
colors=c("Data"="blue","Model"="red","Implied Truth"="#BF40BF")
plot.munic<-ggplot()+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=dat.quants,aes(x=years,y=median,color="Model"))+
  geom_line(data=true.dat.quants,aes(x=years,y=median,color="Implied Truth"),linetype="dashed")+
  geom_line(data=overall_data,aes(x=years,y=median,color="Data"))+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=15))+
  geom_ribbon(data=overall_data,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.4, fill="blue")+
  geom_ribbon(data=dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.2, fill="red")+
  geom_ribbon(data=true.dat.quants,aes(ymin=lowerCI,ymax=upperCI,x=years),alpha=0.08, fill="#BF40BF")+
  scale_color_manual(values = colors,breaks=c("Data","Model","Implied Truth"),limits=c("Data","Model","Implied Truth"))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  xlab("Evolutionary Time (years)")+
  theme(legend.position = "none")+
  ylab("Distance (km)")#+
# plot.munic

###Figure 1E/F if inset is in supplement
library(directlabels)
load("./Figures/Figure1/melted.vaxamr.props.NVT.RData")
load("./Figures/Figure1/melted.vaxamr.props.amr.RData")
# melted.vaxamr.props.NVT<-subset(melted.vaxamr.props, melted.vaxamr.props$type=="NVT")
nvt.props<- ggplot(melted.vaxamr.props.NVT)+
  geom_line(aes(x=year,y=proportion,color="black"))+
  theme_classic()+
  geom_dl(aes(label = type,x=year,y=proportion), method = list(dl.combine("last.points")), cex = 0.8) +
  scale_color_manual(values=c("red"),limits=c("NVT"))+
  theme(legend.position = "none",
        axis.text.x = element_blank(),axis.title = element_text(size=15),axis.title.x = element_blank(),axis.text.y = element_text(size=15))+
  # theme(aspect.ratio=2/6)+
  geom_vline(xintercept = 2009,linetype="dashed",color="black")+
  geom_vline(xintercept = 2011,linetype="dashed",color="black")+
  xlab("Years")+
  xlim(2000,2016)+
  ylab("Proportions")+
  ylim(0,1)
# nvt.props

amr.props<-ggplot(melted.vaxamr.props.amr)+
  theme_classic()+
  geom_line(aes(x=year,y=proportion,group=type),colour="black")+
  geom_dl(aes(label = type,x=year,y=proportion), method = list(dl.combine("last.points")), cex = 0.8) +
  # scale_color_manual(values=c("black","black","black","black"),limits=c("PEN_R","ERY_R","COT_R","CLI_R"))+
  theme(legend.position = "none",
        axis.text = element_text(size=15),axis.title = element_text(size=15))+
  geom_vline(xintercept = 2009,linetype="dashed",color="black")+
  geom_vline(xintercept = 2011,linetype="dashed",color="black")+
  # theme(aspect.ratio=6/6)+
  xlab("Years")+
  xlim(2000,2016)+
  ylab("Proportions")+
  ylim(0,1)
# amr.props

library(patchwork)

(B1+C1)
  (plot.munic/nvt.props/amr.props)


library(ape)
tree5<-read.tree("./Figures/Figure1/tree_gpsc5.nwk")
tree5$tip.label<-rep(NA,length(tree5$tip.label))
tree1<-read.tree("./Figures/Figure1/tree_gpsc1.nwk")
tree1$tip.label<-rep(NA,length(tree1$tip.label))
plot(tree1)
plot(tree5)

