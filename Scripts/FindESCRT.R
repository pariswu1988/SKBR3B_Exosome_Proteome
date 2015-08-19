#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla
#This script assumes you have already run Analysis.R in this session.

# Import the libraries we will use
library("MSnbase")
library("pRoloc")
library("ggplot2")

base <- "C:/WEF_Data/141208_Exosome_Profile_Analysis_Finalv2/"

ClusterList <- ms2df(ClusterData,c("markers","Description","svm","svm.scores"))
GeneNames <- ClusterList$Description

ESCRT <- read.csv(paste0(base,"ESCRT_Proteins.csv"))
ToMatch <- as.vector(ESCRT$Gene)

ContainsTag <- grepl("_",ToMatch)

for (i in 1:length(ToMatch)) {
  if (ContainsTag[i]){
    ToMatch[i] <- gsub("_.*","", ToMatch[i])
  }
  ToMatch[i] <- paste0("GN=",ToMatch[i])
}

MatchIdx <- grepl(paste(ToMatch,collapse="|"),GeneNames)
GeneMatches <- GeneNames[MatchIdx]
FullMatches <- ClusterList[MatchIdx,]

for (n in 1:length(ToMatch)) {
  Match <- grepl(ToMatch[n], FullMatches$Protein.description)
  FullMatches$ESCRT.Protein[Match] <- as.vector(ESCRT$Protein[n])
  FullMatches$ESCRT.Complex[Match] <- as.vector(ESCRT$ESCRT.complex[n])
}

write.csv(FullMatches, paste0(base,"ESCRT_Matches.csv"))