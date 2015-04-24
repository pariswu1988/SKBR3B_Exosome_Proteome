#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

##############################################################################################
### ProteinQuant.R ###########################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. It is designed to calculate protein ratios # 
# from the peptide ratios that were output from 2 replicate experiments in QuantiMORE. This  #
# script relies on ProteinQuant_Functions.R to run.                                          #
#                                                                                            #
##############################################################################################


# Import the libraries we will use
library("plyr")
library("data.table")
source("ProteinQuant_Functions.R")


###############################################################################################
# Import peptide lists results from Quantimore ################################################
###############################################################################################

# "noC" in the file name indicates that I replace commas with semi-colons using as search and 
# replace in Excel. This was necessary to read the file in correctly.
peptides1 <- read.delim("all_peptide_Rep1_noC.txt", header = T)

# Change some columns to more usable names
peptides1 <- rename(peptides1, c("Ratio_130.1411.129.1378_Weighted.Average"="Ratio_130.129",
                                  "Ratio_131.1382.129.1378_Weighted.Average"="Ratio_131.129",
                                  "Ratio_129.1378.129.1378_Weighted.Average"="Ratio_129.129",
                                  "Protien.description"="Description"))

peptides2 <- read.delim("all_peptide_Rep2_noC.txt", header = T)

peptides2 <- rename(peptides2, c("Ratio_130.1411.129.1378_Weighted.Average"="Ratio_130.129",
                                 "Ratio_131.1382.129.1378_Weighted.Average"="Ratio_131.129",
                                 "Ratio_129.1378.129.1378_Weighted.Average"="Ratio_129.129",
                                 "Protien.description"="Description"))

peptides1$Protein.Id <- as.character(peptides1$Protein.Id)
peptides2$Protein.Id <- as.character(peptides2$Protein.Id)
peptides1$Net.peptide <- as.character(peptides1$Net.peptide)
peptides2$Net.peptide <- as.character(peptides2$Net.peptide)

###############################################################################################
# Perform weighted average of Peptides ########################################################
###############################################################################################

#### Combine instances of the same peptide sequence ####
unique_pep1 <- unique(peptides1$Net.peptide) # the sequences of every unique peptide
unique_pep2 <- unique(peptides2$Net.peptide)

# choose columns to keep
ddrop <- c("File.name","nonlabel.peptide","Protein.Id","Description","Net.peptide",
           "Identified.peptide","Score.1","ID.Scan.number","Ratio_130.129","Ratio_131.129","Number")

#Changing the format to data.table allows my terribly inefficient code to run a little faster
peptides1_real <- data.table()
peptides2_real <- data.table()

# These for loops actually do the combining of peptides so that there is only one entry per
# per unique peptide sequence per replicate. Anything variables named pb are progress bars I used to
# make sure the loops were still running.
pb1 <- txtProgressBar(min = 0, max = length(unique_pep1), style = 3)
for (i in 1:length(unique_pep1)){
  current <- unique_pep1[i]
  current <- peptides1[peptides1$Net.peptide %in% current,]
  quant130 <- sum(current$Ratio_130.129*current$Number)/sum(current$Number)
  quant131 <- sum(current$Ratio_131.129*current$Number)/sum(current$Number)
  entry <- current[1,]
  entry$Ratio_130.129[1] <- quant130
  entry$Ratio_131.129[1] <- quant131
  entry$Score.1[1]<- sum(current$Score.1)/nrow(current)
  entry$Number <- sum(current$Number)
  entry <- entry[,names(entry) %in% ddrop]
  entry$Replicate[1] <- 1
  peptides1_real <- rbind(peptides1_real,entry)
  setTxtProgressBar(pb1, i)
}

pb2 <- txtProgressBar(min = 0, max = length(unique_pep2), style = 3)
for (i in 1:length(unique_pep2)){
  current <- unique_pep2[i]
  current <- peptides2[peptides2$Net.peptide %in% current,]
  quant130 <- sum(current$Ratio_130.129*current$Number)/sum(current$Number)
  quant131 <- sum(current$Ratio_131.129*current$Number)/sum(current$Number)
  entry <- current[1,]
  entry$Ratio_130.129[1] <- quant130
  entry$Ratio_131.129[1] <- quant131
  entry$Score.1[1]<- sum(current$Score.1)/nrow(current)
  entry$Number <- sum(current$Number)
  entry <- entry[,names(entry) %in% ddrop]
  entry$Replicate <- 2
  peptides2_real <- rbind(peptides2_real,entry)
  setTxtProgressBar(pb2, i)
}

# Once again, we want character vectors for matching instead of factor vectors.
peptides1_real$Protein.Id <- as.character(peptides1_real$Protein.Id)
peptides2_real$Protein.Id <- as.character(peptides2_real$Protein.Id)

# Filter out DECOYs and calculate peptide and protein level FDR. All decoy's are name with
# DECOY at the front of their protein ID.
decoyidx1 <- grepl("decoy",peptides1_real$Protein.Id, ignore.case = T)
decoyidx2 <- grepl("decoy",peptides2_real$Protein.Id, ignore.case = T)

decoy1 <- peptides1_real[decoyidx1,] # Makes dataframes of decoy peptides
decoy2 <- peptides2_real[decoyidx2,]


pep_fdr1 <- nrow(decoy1)/nrow(peptides1_real) # calculate peptide level FDRs (decoys/total)
pep_fdr2 <- nrow(decoy1)/nrow(peptides2_real)
# These peptide level FDRs are very close to what was calculated by PeptideProphet.


###########################################################################################
####### Building Proteins from Protein Prophet Results ####################################
###########################################################################################

# Parses in the ProteinProphet results saved in a csv file directly from ProteinProphet.
if (!exists("ProtProph_Combined")){
  ProtProph_Combined <- Parse_ProtProph("ProteinProphet_Combined.csv")
}

# Creates list of dataframes containging peptides. Each dataframe contains the peptides
# supporting a single protein according to the ProteinProphet results. These are indexed
# in the same order as the proteins from ProtProph_Combined.
Assigned_Pep_Rep1 <- Assign_Peptides(peptides1_real,ProtProph_Combined)
Assigned_Pep_Rep2 <- Assign_Peptides(peptides2_real,ProtProph_Combined)

# Merges the peptides from both replicates into a List of single dataframes containing
# the peptides quantified in each replicate.
Merged_Peptides <- Merge_Peptides(ProtProph_Combined,Assigned_Pep_Rep1,Assigned_Pep_Rep2)

# Calculates the number of unique peptides used to quantitate the protein
ProtProph_Combined <- Num_Peptides(ProtProph_Combined,Assigned_Pep_Rep1,
                                   Assigned_Pep_Rep2)

# Protein quantitation is performed by replicate.
Proteins_Rep1 <- Prot_Quant(ProtProph_Combined,Assigned_Pep_Rep1)
Proteins_Rep2 <- Prot_Quant(ProtProph_Combined,Assigned_Pep_Rep2)

# Merges the protein lists for both replicates, discarding any proteins observed in
# only 1 replicate.
Proteins_Combined <- merge(Proteins_Rep1,Proteins_Rep2,
                           by = c("Entry_Num","Protein","Probability",
                                  "Coverage","Independent_Spectra",
                                  "Percent_Spectrum_ID","Protein_Group",
                                  "Unique_Peptides"),
                           suffixes = c("_Rep1","_Rep2"))

# Filters the combined protein list to remove decoy proteins and proteins with 
# 0 probability.
Proteins_Combined <- Proteins_Combined[!grepl("^DECOY.*",Proteins_Combined$Protein),]
Proteins_Combined <- Proteins_Combined[Proteins_Combined$Probability != 0,]

# Reads the UniProt identifiers from the database we used and matches them to the 
# protein description.
Key <- ID_Key("UniProtKB-human-2014_DECOY.fasta")
Proteins_Combined <- merge(Proteins_Combined, Key, by = "Protein")

# Calculates the average for each protein between the 2 replicates
Proteins_Combined$Ratio_130.129 <- (Proteins_Combined$Ratio_130.129_Rep1 + 
                                      Proteins_Combined$Ratio_130.129_Rep2)/2
Proteins_Combined$Ratio_131.129 <- (Proteins_Combined$Ratio_131.129_Rep1 + 
                                      Proteins_Combined$Ratio_131.129_Rep2)/2

# Calculates the standard deviation between the protein ratios of the 2 replicates
Proteins_Combined$SD_130.129 <- Col_SD(Proteins_Combined$Ratio_130.129_Rep1,
                                       Proteins_Combined$Ratio_130.129_Rep2)

Proteins_Combined$SD_131.129 <- Col_SD(Proteins_Combined$Ratio_131.129_Rep1,
                                       Proteins_Combined$Ratio_131.129_Rep2)

# again replaces commas with semi-colons to avoid parsing isssues.
Proteins_Combined$Description <- gsub(",",";",Proteins_Combined$Description)

# Matches the identified peptides witha  protein ID and description then saves
# it as a raw peptide list.
Peptides <- Peptide_List(Proteins_Combined,Merged_Peptides)
write.csv(Peptides,"peptide_list.csv", row.names=F)


# Write the Protein results to a CSV.
write.csv(Proteins_Combined,"Proteins.csv", row.names=F)

