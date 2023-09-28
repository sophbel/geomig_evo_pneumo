###############
###Figure 2####
##############
setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
disease=FALSE
sub=TRUE
load("./data/gps_metadata/GPS_SA.RData")###This loads all 6910 SA metadata
if(sub==FALSE){
if(disease==TRUE){
  load("./RelativeRisk/files/mat_lineage.disease.noSub.Rdata")###This loads all 6910 SA metadata
  mat_lineage.noSub<-mat_lineage
  load("./RelativeRisk/matplot_specificRange.disease.noSub.RData")###This loads all 6910 SA metadata
  matplot_specificRange.noSub<-matplot_specificRange
  load("./RelativeRisk/matplot_continuous.disease.noSub.RData")###This loads all 6910 SA metadata
  matplot_continuous.noSub<-matplot_continuous

}else
{
  load("./RelativeRisk/files/mat_lineage.noSub.Rdata")###This loads all 6910 SA metadata
  mat_lineage.noSub<-mat_lineage
  load("./RelativeRisk/files/matplot_specificRange.noSub.RData")###This loads all 6910 SA metadata
  matplot_specificRange.noSub<-matplot_specificRange
  load("./RelativeRisk/files/matplot_continuous.noSub.RData")###This loads all 6910 SA metadata
  matplot_continuous.noSub<-matplot_continuous

}
}else{
  if(disease==TRUE){
  load("./RelativeRisk/files/mat_lineage.disease.Sub.Rdata")###This loads all 6910 SA metadata
  mat_lineage.noSub<-mat_lineage
  load("./RelativeRisk/files/matplot_specificRange.disease.Sub.RData")###This loads all 6910 SA metadata
  matplot_specificRange.noSub<-matplot_specificRange
  load("./RelativeRisk/files/matplot_continuous.disease.Sub.RData")###This loads all 6910 SA metadata
  matplot_continuous.noSub<-matplot_continuous
  # write.table(matplot_continuous.noSub,"continuousDisease_RRs.csv",sep=",",row.names = FALSE,col.names =TRUE)
  
}else
{
  load("./RelativeRisk/files/mat_lineage.Sub.Rdata")###This loads all 6910 SA metadata
  mat_lineage.noSub<-mat_lineage
  load("./RelativeRisk/files/matplot_specificRange.Sub.RData")###This loads all 6910 SA metadata
  matplot_specificRange.noSub<-matplot_specificRange
  # load("matplot_continuous.Sub.RData")###This loads all 6910 SA metadata
  # matplot_continuous.noSub<-matplot_continuous
  # write.table(matplot_continuous.noSub,"continuous_RRs.csv",sep=",",row.names = FALSE,col.names =TRUE)
  
}
}




library(ggplot2)
mat_lineage.noSub$time_range[which(mat_lineage.noSub$time_range=="0 - 1")]<-"Same Lineage"
mat.tmp <- rbind(matplot_specificRange.noSub[,c("distance_range","RR","lowerCI","upperCI","time_range")],
                 mat_lineage.noSub[,c("distance_range","RR","lowerCI","upperCI","time_range")])
mat.tmp$lowerCI[mat.tmp$lowerCI <= 0.1 ] <- 0.1
mat.tmp$upperCI[mat.tmp$upperCI <= 0.1]<- 0.1
mat.tmp$RR[mat.tmp$RR <= 0.1]<- 0.1
mat.tmp$lowerCI[mat.tmp$lowerCI >= 10 ] <- 10
mat.tmp$upperCI[mat.tmp$upperCI >= 10 ] <- 10
mat.tmp$RR[mat.tmp$RR >= 10 ] <- 10

mat.tmp$time_range_f = factor(mat.tmp$time_range, levels=c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"), labels = c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"))

mat.tmp$time_range_f = factor(mat.tmp$time_range, levels=c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"), labels = c("Same Lineage", "0 - 5","5 - 10","10 - 20" ,"20 - 200"))
mat.tmp$distance_range_f = factor(mat.tmp$distance_range, levels=c("Within Province", "<500","500-1000","Distant Pairs" ,"Other Africa","Outside Africa"))
# pdf(file="/Users/sb62/Documents/Migration/SA_Migration_110422/Submission/Figures/Figure_panel2/RRoverDistTime.PDF.pdf",width=3.4,height=10)
mat.tmp$upperCI[which(mat.tmp$upperCI>5)]<-5
lineplot <- ggplot(data = mat.tmp, aes( x = distance_range_f, y = RR , group = time_range_f)) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), alpha = .9, color = "black", width = 0.5 ) +
  theme(axis.text.x = element_text(angle = 90, size=20),
        axis.text.y = element_text(size=18), axis.title.y  = element_text(size=25),
        axis.title.x=element_text(size=18), title =element_text(size=18),
        panel.grid.major= element_blank(), panel.grid.minor= element_blank(),axis.line = element_line(colour = "black")) +
  scale_x_discrete(name ="",
                   
                   limits=c("Within Province", "<500","500-1000","Distant Pairs" ,"Other Africa","Outside Africa"))+
  # theme(aspect.ratio=20/30) +
  
  geom_point(aes(x=4, y=1),shape=25,fill="white", size=1,alpha=1)+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_y_continuous( trans = "log10", breaks = c(0.09,1,5),labels=c("<0.1","1.0","5.0"),limits=c(0.09,5)) +
  ylab("Risk Ratio")

lin1_mal_cont <- lineplot + geom_rect(fill = '#00539CFF', xmin = 0, xmax = 1.5, ymin =-5, ymax = 7, alpha =0.05) +
  geom_rect(fill = '#ED2B33FF', xmin = 1.5, xmax = 4.5, ymin =-5, ymax = 7, alpha =0.05) +   
  geom_rect(fill = '#97BC62FF', xmin = 4.5, xmax = 5.5, ymin =-5, ymax = 7, alpha =0.05) +
  geom_rect(fill = 'purple', xmin = 5.5, xmax = 6.5, ymin =-5, ymax = 7, alpha =0.05) #+
pmanycont <- lin1_mal_cont + facet_wrap( ~ time_range_f, nrow =5) +
  theme(strip.text.x=element_blank())
pmanycont
# dev.off()
# ggsave(pmanycont,file="./RelativeRisk/plots/RR.overdistTime.pdf",width=3.4,height=10)
##### Plot Continuous
mat_distancerange <- matplot_continuous
mat_distancerange$lowerCI[mat_distancerange$lowerCI <= 0.1 ] <- 0.1
mat_distancerange$upperCI[mat_distancerange$upperCI <= 0.1]<- 0.1
mat_distancerange$RR[mat_distancerange$RR <= 0.1]<- 0.1

mat_sub_genomic <- subset(mat_distancerange, mat_distancerange$distance_range == "Within Province" | mat_distancerange$distance_range == "500-1000" | mat_distancerange$distance_range == "Other Africa"| mat_distancerange$distance_range == "Outside Africa")
mat_sub_genomic$distance_range_f = factor(mat_sub_genomic$distance_range, levels= c("Within Province","500-1000","Other Africa","Outside Africa"))
mat_sub_genomic.tmp <- subset(mat_distancerange, mat_distancerange$distance_range=="Within Province")
plotwithin <- ggplot(data = mat_sub_genomic, aes( x = medMRCA, y = RR, group = distance_range_f ,color = distance_range_f)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI), alpha = 0.1, color = NA, fill = "#00539CFF" ) +
  theme(axis.text.x = element_text(angle = 90, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        title = element_text(size=20)) +
  # ggtitle("B")+
  # scale_x_discrete(name ="Distance between isolates (km)",
  # limits=c(matbind$distance_range[1:4]))+
  # scale_y_log10(labels = function(x) format(x, scientific = FALSE))+
  scale_y_continuous( trans = "log10", breaks = c(0.10,1,6),labels=c("<0.1","1.00",">6.00"),limits=c(0.1,6)) +
  # scale_y_continuous( trans = "log10", breaks = c(0.2,1,6),labels=c("<0.2","1.00",">6.00"),limits=c(0.19,6.2)) +
  
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  ylab("Risk Ratio") + 
  xlab("tMRCA (years)")+
  theme(aspect.ratio=14/15) +
  scale_color_manual(limits = c("Within Province","500-1000","Other Africa","Outside Africa"), values = c('#00539CFF','#ED2B33FF','#97BC62FF','purple')) + 
  theme(panel.grid.major =  element_line(colour="lightgrey", size=0.1), panel.grid.minor =  element_line(colour="lightgrey", size=0.1),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  theme(axis.text.y = element_text(size= 15) ) +
  xlim(0,75)
homog <- plotwithin + facet_wrap( ~ distance_range_f, nrow =1) +
  theme(strip.text.x = element_blank())


