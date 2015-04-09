#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

# Import the libraries we will use
library("MSnbase")
library("pRoloc")
library("ggplot2")
library("scales")
library("reshape2")

base <- "C:/WEF_Data/141208_Exosome_Profile_Analysis_Finalv2/"

######################################################################################################
### This script is meant to be used after 3Rep_Analysis.R & 3Rep_Analysis_NoMiss.R ###################
######################################################################################################
MarkerIdx1 <- PlotData$markers != "unknown"
Markers <- PlotData[MarkerIdx1,]
Exo <- Markers[Markers$markers == "Exosome",]
Non <- Markers[Markers$markers == "Non-Exosome",]

#Plots Basic Scatter Plots
BP <- qplot(Ratio_130.129,Ratio_131.129,data = PlotData, 
            xlab = "log2 Enrichement by 100K Spin (TMT 130/129)",
            ylab = "log2 Enrichment by Optiprep (TMT 131/129)")

BP <- BP + theme(legend.background = element_rect(colour = "black"),
                 panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 axis.line = element_line(),
                 text=element_text(size=18))

BP <- BP + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))


#Plots Basic Scatter Plot with Markers
BPM <- ggplot(PlotData, aes(x=Ratio_130.129,y=Ratio_131.129)) + geom_point(color = "#333333") + theme_bw()

BPM <- BPM + geom_point(data = Markers, aes(x=Ratio_130.129,y=Ratio_131.129, colour=markers), size=4)

BPM <- BPM + theme(legend.background = element_rect(colour = "black"),
                   panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                   axis.line = element_line(),
                   text=element_text(size=18),
                   legend.key = element_blank()) +
      xlab("log2 Enrichement by 100K Spin (TMT 130/129)") +
      ylab("log2 Enrichment by Optiprep (TMT 131/129)")

BPM <- BPM + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))
BPM <- BPM + theme(legend.justification=c(0,0), legend.position=c(0,0),
                   axis.text.x = element_text(size=14, colour=rgb(0,0,0)),
                   axis.text.y = element_text(size=14, colour=rgb(0,0,0)))
BPM <- BPM + scale_colour_discrete(name = "Markers")



#Plots Scatter with Clustering
CP <- qplot(Ratio_130.129,Ratio_131.129,data = PlotData, 
            color = svm, size = svm.scores, alpha = I(0.7),
            xlab = "log2 Enrichement by 100K Spin (TMT 130/129)",
            ylab = "log2 Enrichment by Optiprep (TMT 131/129)") + theme_bw()


CP <- CP + theme(legend.background = element_rect(colour = "black"),
                 panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 axis.line = element_line(),
                 text=element_text(size=18),
                 legend.key = element_blank())

CP <- CP + scale_colour_discrete(name = "Prediction") +
  scale_size_continuous(name = "SVM\nProbability", range = c(1,4))

CP <- CP + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))

CPM <- CP + geom_point(data = Markers, aes(x=Ratio_130.129,y=Ratio_131.129, fill=markers),
                       color="black",shape=21) + guides(fill=FALSE)

#CPM


CP2 <- ggplot(data = PlotData, aes(x = Ratio_130.129, y = Ratio_131.129)) + 
  geom_point(aes(fill=svm, shape = svm, size = svm.scores), alpha = I(0.5))+ theme_bw()


CP2 <- CP2 + theme(legend.background = element_rect(colour = "black"),
                panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                axis.line = element_line(),
                text=element_text(size=18),
                legend.key = element_blank()) +
              xlab("log2 Enrichement by 100K Spin (TMT 130/129)") +
              ylab ("log2 Enrichment by Optiprep (TMT 131/129)")

CP2 <- CP2 + scale_size_continuous(name = "SVM\nProbability", range = c(1,4)) + 
  scale_fill_discrete(name="Prediction") +
  scale_shape_manual(name = "Prediction",values=c(21,23,21,23))


CPM2 <- CP2 + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))

#CPM2 <- CPM2 + geom_point(data = Markers, aes(x=Ratio_130.129,y=Ratio_131.129, shape=svm),
#                       color="white", color = "black",size = 4) + guides(fill=FALSE) 

CPM2

Hist <- ggplot(PlotData) +
  geom_density(aes(x=Ratio_130.129), fill="#0099FF",alpha=0.3,color="#0099FF",binwidth=0.2) +
  geom_density(aes(x=Ratio_131.129),fill="#FF9900",alpha=0.3,color="#FF9900",binwidth=0.2) +
  geom_vline(aes(xintercept=0)) +
  xlab("Enrichment (Log2)") +
  ylab("Fraction of Peptides") + scale_y_continuous(expand=c(0,0))
