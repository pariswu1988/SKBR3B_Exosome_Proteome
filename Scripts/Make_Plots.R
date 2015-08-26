##############################################################################################
### Make_Plots.R #############################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. It is designed to calculate protein ratios # 
# from the peptide ratios that were output from 2 replicate experiments in QuantiMORE.       #
# exososomes using the TMT 129, 130 and 131 tags. This script makes the plots.               #
#                                                                                            #
##############################################################################################



library("plyr")
library("ggplot2")
library("grid")
source("Scripts/Plot_Functions.R")

dir.create("Figures")

data <- read.csv("Temp/ClusterResults.csv")
data <- rename(data,c("X" = "Protein.Id"))


Cluster <- cluster_plot_marked(data)
Cluster
ggsave(file = "Figures/ClusterPlot.pdf",width=8.5,height=8.10,units = "cm", useDingbats=F)
ggsave(file = "Figures/ClusterPlot.tiff",width=8.5,height=8.10,units = "cm")

#################################################################
cover_fill = c("#FF6666","#66FFFF","#FF6666","#66FFFF")
cover_color = c("black","black","black","black")

Cover <- cluster_plot_marked(data) +
  scale_fill_manual(values = cover_fill, name = "Classification") +
  scale_color_manual(values = cover_color, name = "Classification") +
  scale_shape_manual(values = c(21,23,21,23), name = "Classification") +
  scale_alpha_manual(values = c(0.7,0.7,1,1), name = "Classification") +
  theme(legend.position = "none",axis.title.x = element_text(size=6),axis.title.y = element_text(size=6)) +
  xlab(expression("log"[2]~"Enrichment in 100K x g Ultracentrifugation")) +
  ylab(expression("log"[2]~"Enrichment in Density Gradient"))

Cover
ggsave(file = "Figures/CoverPlot.tiff",width=9, height = 5, units = "cm")

#################################################################

Scatter_M <- scatter_plot_marked(data)
Scatter_M
ggsave(file = "Figures/ScatterPlot_M.pdf",width=8.5,height=6.57,units = "cm",useDingbats=F)
ggsave(file = "Figures/ScatterPlot_M.tiff",width=8.5,height=6.57,units = "cm")

#################################################################

Scatter_M_Helper <- scatter_plot_marked(data)
markers <- read.csv("Data/markers.csv", as.is = T)
data2 <- data
data2$Protein.Id <- levels(data2$Protein.Id)
markers$Regex <- sapply(markers$Protein, function(x) gsub("\\|","\\\\|",x))
markers$R131 <- sapply(markers$Regex, function(x) data$Log2_Ratio_131.129[grepl(x,data2$Protein.Id)])
markers$R130 <- sapply(markers$Regex, function(x) data$Log2_Ratio_130.129[grepl(x,data2$Protein.Id)])
markers$ID <- gsub("^.*\\|","",markers$Protein)

Scatter_M_Helper <- Scatter_M_Helper +  
  annotate("text", x=markers$R130, y=markers$R131, label=markers$ID, colour="red")

Scatter_M_Helper
  
ggsave(file="Figures/ScatterPlot_Help.pdf",width=85,height=65.7,units="cm",useDingbats=F)

##################################################################

Scatter <- scatter_plot(data)
Scatter
ggsave(file = "Figures/ScatterPlot.pdf",width=8.5,height=6.57,units = "cm",useDingbats=F)
ggsave(file = "Figures/ScatterPlot.tiff",width=8.5,height=6.57,units = "cm")

###################################################################

PM_Markers <- cluster_plot_om(data)
PM_Markers
ggsave(file = "Figures/PM_MarkerPlot.pdf",width=8.5,height=8.07,units = "cm",useDingbats=F)
ggsave(file = "Figures/PM_MarkerPlot.tiff",width=8.5,height=8.07,units = "cm")

###################################################################
PM_csv <- map_markers(data)
PM_csv <- PM_csv[!is.na(PM_csv$pRolocmarkers),]
PM_csv2 <- PM_csv[,c(2,3,4,5,6,7,9)]
PM_csv2 <- rename(PM_csv, c("Ratio_130.129" = "log2_Ratio_130.129",
                           "Ratio_131.129" = "log2_Ratio_131.129"))
write.csv(PM_csv2,"Temp/PM_markers.csv",row.names=F)

###################################################################

PM_Markers_Helper <- cluster_plot_om(data)
PM_Markers_Helper <- PM_Markers_Helper + 
  annotate("text",x=PM_csv$Log2_Ratio_130.129,y=PM_csv$Log2_Ratio_131.129,label=PM_csv$Gene_Name, color="red")
PM_Markers_Helper

ggsave(file="Figures/PM_Marker_Help.pdf",width=85, height=80.7, units="cm",useDingbats=F)


###################################################################
data$Log2_Ratio_131.130 <- log2(data$Ratio_131.129/data$Ratio_130.129)

revPlot <- cluster_plot_rev(data)
revPlot

ggsave(file = "Figures/Cluster_rev.tiff",width=8.5,height=8.10,units = "cm")


## Validation Venn Diagram #######################################
source("Scripts/Validation.R")
