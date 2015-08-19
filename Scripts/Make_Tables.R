##############################################################################################
### Make_Tables.R ############################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. This script is used to create tables from  #
# the cluster results. Used after ClusterAnalysis.R and Make_Plots.R.                        #
#                                                                                            #
##############################################################################################


library(plyr)
source("Scripts/ProteinQuant_Functions.R")

dir.create("Tables")

Proteins <- read.csv("Temp/ClusterResults.csv")

#Proteins <- rename(Proteins, replace = c("X" = "Protein_ID"))
names(Proteins)[1] <- "Protein_ID"

# Fix Accession numbers:
accession <- gsub("^.*?\\|","",Proteins$Protein_ID)
accession <- gsub("\\|.*$","",accession)

Proteins$accession <- accession

# Rounding Ratios
Proteins$Ratio_130.129 <- format(round(Proteins$Ratio_130.129,3),nsmall=3)
Proteins$Ratio_131.129 <- format(round(Proteins$Ratio_131.129,3),nsmall=3)
Proteins$SD_130.129 <- format(round(Proteins$SD_130.129,3),nsmall=3)
Proteins$SD_131.129 <- format(round(Proteins$SD_131.129,3),nsmall=3)
Proteins$Log2_Ratio_130.129 <- format(round(Proteins$Log2_Ratio_130.129,3),nsmall=3)
Proteins$Log2_Ratio_131.129 <- format(round(Proteins$Log2_Ratio_131.129,3),nsmall=3)
Proteins$svm.scores <- format(round(Proteins$svm.scores,3),nsmall=3)

# Create Ratios +/- SD
Proteins$Ratio_130_SD <- paste(Proteins$Ratio_130.129,Proteins$SD_130.129, sep = " +/- ")
Proteins$Ratio_131_SD <- paste(Proteins$Ratio_131.129,Proteins$SD_131.129, sep = " +/- ")

Proteins <- Prot_Group(Proteins)

########### Supplmentary Table 1: All Proteins #####################################################
TS1names <- c("accession","Description","Gene_Name","Protein_Group","Ratio_130_SD","Log2_Ratio_130.129",
              "Ratio_131_SD","Log2_Ratio_131.129","svm","svm.scores")
TableS1 <- Proteins[,TS1names]

names(TableS1) <- c("Uniprot Accession","Protein Description","Gene Symbol","Protein Group", 
                    "TMT 130/129 +/- SD","Log2 TMT 130/129","TMT 131/129 +/- SD",
                    "Log2 TMT 131/129","SVM Classification","SVM Probability")

write.csv(TableS1,"Tables/table_S1_AllProteins.csv",row.names=F)

########## Table 1: Exosome and Non-Exosome Markers ################################################
T1names <- c("Gene_Name", "Description","Log2_Ratio_130.129","Log2_Ratio_131.129","Unique_Peptides",
             "markers")
Table1 <- Proteins[Proteins$markers == "Exosome" | Proteins$markers == "Non-Exosome",T1names]

names(Table1) <- c("Gene Name","Description","Log2 TMT 130/129","Log2 TMT 131/129","Unique Peptides",
                   "Marker Class")

write.csv(Table1,"Tables/table_1_ExosomeMarkers.csv", row.names=F)

########## NOT BEING DONE: Top Ten Exosome and Non Exosome #########################################
# ### NOTE: I removed "(Fragment)" from protein description in this case
# T2names <- c("Gene_Name","Description","Log2_Ratio_130.129","Log2_Ratio_131.129",
#              "svm.scores","svm")
# #NoFrag <- baseTable
# #NoFrag$Description <- gsub(" \\(Fragment\\)$","",NoFrag$Description) 
# Exo <- Proteins[Proteins$svm == "Exosome" & Proteins$markers == "unknown",]
# Exo <- arrange(Exo, desc(svm.scores))
# Exo <- Exo[1:10,T2names]
# write.csv(Exo,"Tables/table_2_Exo.csv", row.names=F)
# 
# NonExo <- Proteins[Proteins$svm == "Non-Exosome" & Proteins$markers == "unknown",]
# NonExo <- arrange(NonExo, desc(svm.scores))
# NonExo <- NonExo[1:10,T2names]
# write.csv(NonExo,"Tables/table_2_NonExo.csv", row.names=F)
# 
# Table2 <- rbind(Exo,NonExo)
# write.csv(Table2,"Tables/table_2_both.csv", row.names=F)

######### Table 3: PM marker table ###################################################################
PM <- read.csv("Temp/PM_Markers.csv")
PM <- PM$Protein.Id
TS2names <- c("svm","Gene_Name","Description","Log2_Ratio_130.129","Log2_Ratio_131.129","svm.scores")
TableS2 <- Proteins[Proteins$Protein_ID %in% PM, TS2names]
TableS2 <- arrange(TableS2,svm)

names(TableS2) <- c("SVM Classification", "Gene Name", "Description","Log2 TMT 130/129", "Log2 TMT 131/129",
                    "SVM Probability")

write.csv(TableS2,"Tables/table_S2_PMMarkers.csv", row.names=F)



#####################################################################################################
######### Supplementary Information #################################################################
#####################################################################################################

##### SI list 1: Identified Proteins ################################################################
SI1names <- c("Protein_Group","accession","Gene_Name","Description","Probability","svm","svm.scores",
              "Ratio_130_SD","Log2_Ratio_130.129","Ratio_131_SD","Log2_Ratio_131.129", "Unique_Peptides",
              "Coverage")
SI1 <- Proteins[,SI1names]
names(SI1) <- c("Protein Group", "Uniprot Accession", "Gene Symbol", "Protein Description",
                "ProteinProphet Probability", "SVM Classification", "SVM Probability","TMT 130/129 +/- SD",
                "Log2 TMT 130/129", "TMT 131/129 +/- SD", "Log2 TMT 131/129", "Unique Peptides",
                "Coverage (%)")
write.csv(SI1,"Tables/SI-1_Proteins.csv",row.names=F)

##### SI list 2: Peptides ###########################################################################
Peptides <- read.csv("Temp/peptide_list.csv", as.is = T)
SI2names <- c("Assigned_Protein","Net.peptide","Ratio_130.129","Ratio_131.129","Number","Replicate")
SI2 <- Peptides[,SI2names]
SI2 <- SI2[SI2$Assigned_Protein %in% Proteins$Protein_ID,]

Match_ID <- function(x,field) {
  b <- Proteins[x == Proteins$Protein_ID,field]
  return(b)
}

#SI2$Desc <- sapply(SI2$Assigned_Protein, Match_ID,"Description")
SI2$Acc <- sapply(SI2$Assigned_Protein, Match_ID,"accession")
SI2$GN <- sapply(SI2$Assigned_Protein, Match_ID,"Gene_Name")

SI2order <- c("Acc","GN","Net.peptide","Ratio_130.129","Ratio_131.129","Number","Replicate")
SI2 <- SI2[,SI2order]

names(SI2) <- c("Uniprot Accession","Gene Name","Peptide Sequence","TMT 130/129","TMT 131/129",
                "PSMs","Replicate")

write.csv(SI2,"Tables/SI-2_Peptides.csv",row.names=F)



##### SI list 3: KRONA Classifications ##############################################################
# Generate gene lists to do Krona Classifications                                                   #
#####################################################################################################

GN_Exosome <- Proteins$Gene_Name[Proteins$svm == "Exosome"]
GN_Non_Exosome <- Proteins$Gene_Name[Proteins$svm == "Non-Exosome"]
write.table(GN_Non_Exosome,"Temp/GN_Non.txt",row.names=F,quote=F,col.names=F)
write.table(GN_Exosome,"Temp/GN_Exo.txt",row.names=F,quote=F,col.names=F)

#### Import the Results from DAVID ##################################################################
DAVID_Exo <- read.delim("Data/DAVID_Exo.txt", as.is=T)
DAVID_Exo$GO <- gsub("~.*$","",DAVID_Exo$Term)
DAVID_Non <- read.delim("Data/DAVID_Non_Exo.txt", as.is=T)
DAVID_Non$GO <- gsub("~.*$","",DAVID_Non$Term)
GoTerms <- read.csv("Data/GoTerms.csv",as.is=T)
Exo <- sapply(GoTerms$GO.ID,function(x) DAVID_Exo$Count[DAVID_Exo$GO == x])
NonExo <- sapply(GoTerms$GO.ID,function(x) DAVID_Non$Count[DAVID_Non$GO == x])

for(i in 1:nrow(GoTerms)){
  if(length(Exo[[i]]) == 0 ){ Exo[[i]] <- 0 }
  if(length(NonExo[[i]]) == 0 ){ NonExo[[i]] <- 0}
  
}

GoTerms$Exo <- unlist(Exo)
GoTerms$NonExo <- unlist(NonExo)

write.csv(GoTerms, "Tables/Krona_GO_Counts.csv",row.names=F)


for(i in 1:nrow(GoTerms)){
  a <- DAVID_Exo$Gene[DAVID_Exo$GO == GoTerms$GO.ID[i]]
  b <- DAVID_Non$Gene[DAVID_Non$GO == GoTerms$GO.ID[i]]
  if(length(a)==0){a <- ""}
  if(length(b)==0){b <- ""}
  GoTerms$ExoProt[i] <- a
  GoTerms$NonExoProt[i] <- b
}
Proteins$GO <- rep("",nrow(Proteins))

for(i in 1:nrow(GoTerms)){
  genes <- c(unlist(strsplit(GoTerms$ExoProt[i], ", ", fixed = T)),
             unlist(strsplit(GoTerms$NonExoProt[i],", ", fixed = T)))
  idx <- Proteins$Gene_Name %in% genes
  Term <- paste(GoTerms$GO.ID[i], GoTerms$Name[i], sep = " ")
  Proteins$GO[idx] <- paste(Proteins$GO[idx],Term, sep= "; ")
}

Proteins$GO <- gsub("^; ","", Proteins$GO)

KRONA <- Proteins[,c("accession","Gene_Name","Description","GO","svm")]
names(KRONA) <- c("Uniprot Accession","Gene Name","Description","Gene Ontology Classifications",
                  "SVM Classification")

write.csv(KRONA,"Tables/SI-3_Krona.csv",row.names=F)

##### Panther Classifications #######################################################################
# David should be handling these                                                                    #
#####################################################################################################