##############################################################################################
### Validation.R #############################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. It is designed to calculate protein ratios # 
# from the peptide ratios that were output from 2 replicate experiments in QuantiMORE. This  #
# script compares the validation SEC purification of exosomes with our initial findings      #
#                                                                                            #
##############################################################################################
library("Vennerable")
library("VennDiagram")
library("venneuler")
library("ggplot2")
library("grid")
source("scripts/ProteinQuant_Functions.R")

secProt <- Parse_ProtProph("Data/ProteinProphet_SEC.csv")
secProt <- secProt[!grepl("^DECOY.*",secProt$Protein),]
secProt <- secProt[secProt$Probability != 0,]

Key <- ID_Key("Data/UniProtKB-human-2014_DECOY.fasta")
secProt <- merge(secProt, Key, by = "Protein")

optiProt <- read.csv("Temp/ClusterResults.csv")
optiProt <- optiProt[!is.na(optiProt$svm),]

optiGN <- as.character(optiProt$X)
secGN <- as.character(secProt$Protein)

idx <- optiProt$svm == "Exosome"

exoGN <- optiGN[idx]
nonExoGN <- optiGN[!idx]

exoVal <- sum(exoGN %in% secGN)
nonExoVal <- sum(nonExoGN %in% secGN)



listGN <- list(secGN,exoGN,nonExoGN)
area1 <- length(secGN)
area2 <- length(exoGN)
area3 <- length(nonExoGN)
n12 <-  sum(secGN %in% exoGN)
n23 <- sum(exoGN %in% nonExoGN)
n13 <-  sum(secGN %in% nonExoGN)

data <- data.frame("count" = c(area1,area2,area3), 
                   "x" = c(2.35,3,1.15), 
                   "y" = c(1,1,1),
                   "names" = c("SEC Identifications", 
                               "Exosome Cluster", 
                               "Non-Exosome Cluster"))
                 
plot <- ggplot(data, aes(x=x,y=y,size=count,fill=names)) + 
  geom_point(shape=21, color = "black", alpha = 0.5) + 
  geom_point(shape=1 , color = "black") +
  scale_size_area(max_size= 50) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0),"cm"),
        axis.line = element_blank(),
        panel.border= element_blank(),
        panel.margin = unit(c(0,0,0,0), "cm"),
        panel.background = element_blank()) +
  xlim(c(0,3.5)) +
  ylim(c(0.5,1.5)) +
  annotate("text", x=1.15 , y=0.52, 
           label= paste0("Non-Exosome Cluster\n",area3), size=3, lineheight=0.75) +
  annotate("text", x=2.35, y=1.35, 
           label = paste0("SEC\nIdentifications\n",area1), size=3, lineheight=0.75) +
  annotate("text", x=3, y=0.72, 
           label = paste("Exosome\nCluster\n", area2), size=3, lineheight=0.75) +
  annotate("text", x=1.15, y=1, label = area3 - n13, size=3) +
  annotate("text", x=2.35, y=1, label = area1 - n13 - n12, size=3) +
  annotate("text", x=3.1, y=1, label = area2 - n12, size=3) +
  annotate("text", x=2.78, y=1, label = n12, size=3) +
  annotate("text", x=1.95, y=1, label = n13, size=3)
        

plot

ggsave(file = "Figures/Validation_Venn.tiff",width=8.5,height=6,units = "cm")

write.csv(secProt,"Temp/SEC_proteins.csv",row.names=F)
