#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

##############################################################################################
### Plot_Functions.R ###########################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. This script is designed to create the      #
# plots for the cluster analysis results. Meant to be used after ClusterAnalysis.R           #
#                                                                                            #
##############################################################################################



# Import the libraries we will use
library("ggplot2")
library("scales")
library("reshape2")
library("pRoloc")

my_fill <- c("#FF6666","#66FFFF","black","black")
my_color <- c("black","black","white","white")
xtitle <- expression("Enrichment in 100K Spin (log"[2]~"TMT 130/129)")
ytitle <- expression("Enrichment in Optiprep (log"[2]~"TMT 131/129)")
ytitle2 <- expression("Enrichment in Optiprep (log"[2]~"TMT 131/130)")

scatter_plot <- function(data) {
  plot <- ggplot(data, aes(x=Log2_Ratio_130.129,y=Log2_Ratio_131.129), environment = environment()) + 
    geom_point(color = "#333333", size = 1) + 
    theme_bw() +
    theme(legend.background = element_rect(colour = "black"),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.border = element_rect(colour="black"),
          axis.line = element_line(),
          text=element_text(size=8),
          legend.key = element_blank()) +
    xlab(xtitle) +
    ylab(ytitle) +
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))  
  plot
}

scatter_plot_marked <- function(data){
  markers <- data[data$markers != "unknown",]
  plot <- scatter_plot(data) + 
    geom_point(color = "grey", size=1) +
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
    geom_point(data = markers, 
               aes(x=Log2_Ratio_130.129,y=Log2_Ratio_131.129, fill = markers, shape = markers),
               size = 2, 
               color = "black") +
    scale_fill_manual(values=my_fill, name = "Marker") +
    scale_shape_manual(values = c(21,23), name = "Marker") +
    theme(legend.background = element_blank(),
          legend.justification=c(0.83,0.2), 
          legend.position=c(1,0),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(0,"cm"),
          panel.border = element_rect(colour="black"))
    
  
  plot
}



cluster_plot_marked <- function(data){
  marked <- data
  marked$svm <- as.character(marked$svm)
  marked$svm[marked$markers == "Exosome"] = "Exosome Marker"
  marked$svm[marked$markers == "Non-Exosome"] = "Non-Exosome Marker"
  marked$svm <- as.factor(marked$svm)
  marked$svm <- factor(marked$svm,levels(marked$svm)[c(1,3,2,4)])
  marked <- marked[order(marked$svm),]
  
  plot <- ggplot(marked, aes(x=Log2_Ratio_130.129,y=Log2_Ratio_131.129), environment = environment()) + 
    geom_point(aes(fill = svm, size = svm.scores, shape = svm, alpha=svm, color=svm)) + 
    theme_bw() +
    scale_fill_manual(values = my_fill, name = "Classification") +
    scale_color_manual(values = my_color, name = "Classification") +
    scale_shape_manual(values = c(21,23,21,23), name = "Classification") +
    scale_alpha_manual(values = c(0.7,0.7,1,1), name = "Classification") +
    scale_size_continuous(range=c(0.5,2))+
    theme(legend.background = element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.line = element_line(),
          text=element_text(size=8),
          legend.key = element_blank(),
          legend.justification=c(0,0), 
          legend.position= "bottom",
          legend.title=element_blank(),
          legend.key.size = unit(0,"cm"),
          panel.border = element_rect(colour="black")) +
    xlab(xtitle) +
    ylab(ytitle) +
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
    guides(colour = guide_legend(override.aes = list(size=2)), 
           size = F, fill = guide_legend(ncol=2,byrow=T))
  plot
}

cluster_plot_om <- function(data, marker = "PM") {
  marked <- map_markers(data,marker)
  marked$pRolocmarkers <- as.character(marked$pRolocmarkers)
  marked$svm <- as.character(marked$svm)
  
  for(i in 1:nrow(marked)) {
    if(is.na(marked$pRolocmarkers[i])){
      marked$pRolocmarkers[i] <- marked$svm[i]
    } else {
      marked$svm.scores[i] <- 1
      marked$svm[i] <- paste0(marker, " marker (",marked$svm[i],")")
    }
  }
  marked$pRolocmarkers <- as.factor(marked$pRolocmarkers)
  marked$svm <- as.factor(marked$svm)
  marked <- marked[order(marked$svm),]
  
  plot <- ggplot(marked, aes(x=Log2_Ratio_130.129,y=Log2_Ratio_131.129), environment = environment()) + 
    geom_point(aes(fill = svm, size = svm.scores, shape = svm, alpha=svm, color = svm)) + 
    theme_bw() +
    scale_fill_manual(values = my_fill, name = "Classification") +
    scale_color_manual(values = my_color, name = "Classification") +
    scale_shape_manual(values = c(21,23,21,23), name = "Classification") +
    scale_alpha_manual(values = c(0.7,0.7,1,1), name = "Classification") +
    scale_size_continuous(range=c(0.5,2))+
    theme(legend.background = element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.line = element_line(),
          text=element_text(size=8),
          legend.key = element_blank(),
          legend.justification=c(0,0), 
          legend.position= "bottom",
          legend.title=element_blank(),
          legend.key.size = unit(0,"cm"),
          panel.border = element_rect(colour="black")) +
    xlab(xtitle) +
    ylab(ytitle) +
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
    guides(colour = guide_legend(override.aes = list(size=2)), size = F, 
           size = F, fill = guide_legend(ncol=2,byrow=T))
  plot
}

map_markers <- function(data, mark = ("PM")){
  marker <- mark
  marked <- data
  pm <- data.frame(pRolocmarkers("hsap"))
  pm <- data.frame(row.names(pm),pm)
  pm <- rename(pm, c("pRolocmarkers..hsap.." = "pRolocmarkers"))
  pm <- pm[pm$pRolocmarkers == marker,]
  names(pm)[1] <- "Id"
  
  marked$UPID <- marked$Protein.Id
  marked$UPID <- gsub("^..\\|.*\\|", "", marked$UPID)
  marked <- merge(marked,pm,by.x="UPID", by.y="Id", all.x = T, all.y = F)
  marked
  
}


cluster_plot_rev <- function(data){
  marked <- data
  marked$svm <- as.character(marked$svm)
  marked$svm[marked$markers == "Exosome"] = "Exosome Marker"
  marked$svm[marked$markers == "Non-Exosome"] = "Non-Exosome Marker"
  marked$svm <- as.factor(marked$svm)
  marked$svm <- factor(marked$svm,levels(marked$svm)[c(1,3,2,4)])
  marked <- marked[order(marked$svm),]
  
  plot <- ggplot(marked, aes(x=Log2_Ratio_130.129,y=Log2_Ratio_131.130), environment = environment()) + 
    geom_point(aes(fill = svm, size = svm.scores, shape = svm, alpha=svm, color=svm)) + 
    theme_bw() +
    scale_fill_manual(values = my_fill, name = "Classification") +
    scale_color_manual(values = my_color, name = "Classification") +
    scale_shape_manual(values = c(21,23,21,23), name = "Classification") +
    scale_alpha_manual(values = c(0.7,0.7,1,1), name = "Classification") +
    scale_size_continuous(range=c(0.5,2))+
    theme(legend.background = element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.line = element_line(),
          text=element_text(size=8),
          legend.key = element_blank(),
          legend.justification=c(0,0), 
          legend.position= "bottom",
          legend.title=element_blank(),
          legend.key.size = unit(0,"cm"),
          panel.border = element_rect(colour="black")) +
    xlab(xtitle) +
    ylab(ytitle2) +
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
    guides(colour = guide_legend(override.aes = list(size=2)), 
           size = F, fill = guide_legend(ncol=2,byrow=T))
  plot
}
