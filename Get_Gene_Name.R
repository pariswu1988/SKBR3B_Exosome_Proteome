#! /usr/bin/Rscript --vanilla

base <- "C:/WEF_Data/141208_Exosome_Profile_Analysis_Finalv2/"

ClusterResults <- read.csv(paste0(base,"ClusterResults.csv"))

description <- ClusterResults$Description

GeneName <- gsub("^.*GN=","",description)
GeneName <- gsub(" .*$","",GeneName)

ClusterResults$Gene_Name <- GeneName

ClusterResults_GN <- write.csv(ClusterResults, paste(base,"ClusterResults_GN.csv"))