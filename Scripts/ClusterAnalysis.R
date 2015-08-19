#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

##############################################################################################
### ClusterAnalysis.R ########################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. This script is designed to implement the   #
# pRoloc package to perform SVM classification of our protein list based. This script should #
# be used after ProteinQuant.R.                                                              #
#                                                                                            #
##############################################################################################


# Import the libraries we will use
library("MSnbase")
library("pRoloc")
library("ggplot2")
library("plyr")

##################################################################################################
####### Importing Protein Quantitation Results ###################################################
##################################################################################################

# Importing your Protein Data
proteins <- read.csv("Temp/Proteins.csv")

proteins$Log2_Ratio_130.129 <- proteins$Ratio_130.129
proteins$Log2_Ratio_131.129 <- proteins$Ratio_131.129

write.csv(data.frame(proteins$Protein, 
                     proteins$Log2_Ratio_130.129, 
                     proteins$Log2_Ratio_131.129), 
          "Temp/exprsData.csv", row.names = F,quote = F)

exprsData <- "Temp/exprsData.csv"



write.csv(data.frame(proteins$Protein,
                     proteins$SD_130.129,
                     proteins$SD_131.129,
                     proteins$Description,
                     proteins$Unique_Peptides,
                     proteins$Gene_Name,
                     proteins$Probability,
                     proteins$Coverage,
                     proteins$Entry_Num,
                     proteins$Protein_Group,
                     proteins$Ratio_130.129,
                     proteins$Ratio_131.129),
          "Temp/fData.csv", row.names = F,quote = F)

fData <- "Temp/fData.csv"

##################################################################################################
####### Data Transformation and Clustering #######################################################
##################################################################################################

# Reading your Isoquant Data in for the pRoloc tools
Data <- readMSnSet(exprsFile = exprsData, featureDataFile = fData, sep=",")
Data <- log(Data, 2) # Changed to log2 scale for fold change

#Add localization markers to the data
marked <- addMarkers(Data, "Data/markers.csv")

#Cluster analysis using SVM (for good results time should = 100)
#Seed set for reproducibility.
params <- svmOptimisation(marked, time = 100, xval=5, seed = 1)
ClusterData <- svmClassification(marked,params)

#Save resulting data into a Data frame for easy plotting and manipulation
PlotData <- ms2df(ClusterData,c("proteins.Description","proteins.SD_130.129","proteins.SD_131.129",
                                "markers","svm","svm.scores","proteins.Unique_Peptides",
                                "proteins.Probability","proteins.Entry_Num",
                                "proteins.Coverage","proteins.Gene_Name",
                                "proteins.Protein_Group","proteins.Ratio_130.129",
                                "proteins.Ratio_131.129","proteins.Coverage"))

PlotData <- rename(PlotData, c("proteins.Ratio_130.129"="Ratio_130.129",
                               "proteins.Ratio_131.129"="Ratio_131.129",
                               "proteins.Description"="Description",
                               "proteins.SD_130.129"="SD_130.129",
                               "proteins.SD_131.129"="SD_131.129",
                               "proteins.Unique_Peptides"="Unique_Peptides",
                               "proteins.Probability"="Probability",
                               "proteins.Entry_Num"="Entry_Num",
                               "proteins.Coverage"="Coverage",
                               "proteins.Gene_Name"="Gene_Name",
                               "proteins.Protein_Group"="Protein_Group",
                               "proteins.Log2_Ratio_130.129"="Log2_Ratio_130.129",
                               "proteins.Log2_Ratio_131.129"="Log2_Ratio_131.129",
                               "proteins.Coverage"="Coverage"))

Order <- c("Description","Gene_Name","Probability","svm","svm.scores",
           "Log2_Ratio_130.129","Log2_Ratio_131.129","Ratio_130.129",
           "SD_130.129","Ratio_131.129","SD_131.129","Unique_Peptides",
           "Protein_Group","markers","Coverage")

CSV <- PlotData[,Order] #Reorder columns for spreadsheet
write.csv(CSV, file = "Temp/ClusterResults.csv")
