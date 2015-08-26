#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

##############################################################################################
### ProteinQuant_Functions.R #################################################################
##############################################################################################
#                                                                                            #
# This code was written by Will Fondrie as part of Dr. Austin Yang's lab at the              #
# University of Maryland, Baltimore.                                                         #
#                                                                                            #
# This code was intended for the proteomic analysis of our TMT 6-plex labeled SKBR3b         #
# exososomes using the TMT 129, 130 and 131 tags. The functions in this script are used in   # 
# ProteinQuant.R to perform peptide and protein quantitation.                                #
#                                                                                            #
##############################################################################################

### Reads in ProteinProphet results from the xls (converted to csv) output ###
Parse_ProtProph <- function(filename){  
  
  ProtProph <- read.csv(filename, stringsAsFactors = F)  
  
  print("Parsing ProteinProphet Results")
  b <- data.frame(Entry_Num = character(), Protein = character(), Probability = numeric(),
                  Coverage = numeric(), Independent_Spectra = numeric(),
                  Percent_Spectrum_ID = numeric(), 
                  Peptides = character(), stringsAsFactors = F)
  drop <- "protein.description"
  ProtProph <- ProtProph[,names(ProtProph) != drop]
  
  k <- 1
  pb <- txtProgressBar(min = 0, max = nrow(ProtProph), style = 3)
  
  suppressWarnings(for (i in 1:nrow(ProtProph)) {
    if (!grepl("\\D",ProtProph[i,1]) && ProtProph[i,2]!= ""){
      b[k,] <- as.vector(ProtProph[i,])
      b$Protein_Group[k] <- ProtProph[i,1]
      k <- k + 1
    } else if (grepl("\\D",ProtProph[i,1])) {
      b[k,] <- as.vector(ProtProph[i,])
      b$Protein_Group[k] <- gsub("\\D","",ProtProph[i,1])
      k <- k + 1
    }
    setTxtProgressBar(pb, i)
  }) 
  return(b)
}

### Assigns Peptides to Proteins ###
# Output of this function is the assigned peptides for a protein in the same
# order as the protein list. 
Assign_Peptides <- function(Peptides,Proteins){
  Peptides <- as.data.frame(Peptides)
  
  Pep_Map <- lapply(Proteins$Peptides,strsplit,"\\+")
  
  labels <- c("Net.peptide","Ratio_130.129","Ratio_131.129","Number","Replicate")
  Assigned_Peptides <- lapply(Pep_Map, function(x) Peptides[(Peptides$Net.peptide %in% unlist(x)),labels])
  names(Assigned_Peptides) <- Proteins$Protein
  
  return(Assigned_Peptides)
}

### Merge Peptide Lists ###
Merge_Peptides <- function(ProtList,PepList1,PepList2){
  peptides <- list()
  for(i in 1:nrow(ProtList)){
    peptides[[i]] <- merge(PepList1[[i]],PepList2[[i]],all=T)
  }
  names(peptides) <- ProtList$Protein
  peptides
}

### Protein Quantitation ###
Prot_Quant <- function(Proteins,Assigned_Peptides){
  
  # in this X will be the peptide dataframe, and R refers to an element of 
  # "Ratios" defined below.
  Quant <- function(X, R){
    # each of these performs this type of average:
    # AVERAGE = sum( #PSMs * PeptideRatio ) / Total#PSMs
    avg_130 <- sum(X[,"Number"]*X[,R[1]])/sum(X[,"Number"])
    avg_131 <- sum(X[,"Number"]*X[,R[2]])/sum(X[,"Number"])
    return(c(avg_130,avg_131))
  }
  
  Ratios <- c("Ratio_130.129", "Ratio_131.129")
  Prot_Ratios <- sapply(Assigned_Peptides, Quant, Ratios)
  
  drop <- c("Peptides")
  
  Proteins$Ratio_130.129 <- Prot_Ratios[1,]
  Proteins$Ratio_131.129 <- Prot_Ratios[2,]
  Proteins <- Proteins[is.finite(Proteins[,"Ratio_130.129"]) & is.finite(Proteins[,"Ratio_131.129"]),]
  Proteins <- Proteins[,names(Proteins)!= drop]
  return(Proteins)
}

### Make Protein ID to Description Key ###
ID_Key <- function(Database){
  DB_File <- file(Database)
  DB <- readLines(DB_File)
  close(DB_File)
  
  DB <- DB[grepl("^>", DB) & !grepl("^>DECOY", DB)]
  
  Protein <- gsub("\\s.*$", "",DB)
  Protein <- gsub("^>","",Protein)
  
  Description <- gsub("^.*?\\s","", DB)
  Description <- gsub("OS=H.*$","",Description)
  
  Gene_Name <- gsub("^.*GN=","",DB)
  Gene_Name <- gsub("\\s.*$","",Gene_Name)
  Gene_Name <- gsub("^>.*","",Gene_Name)
  
  Key <- data.frame(Protein = Protein, Description = Description, Gene_Name = Gene_Name)
}

### Standard Deviation Between 2 columns ###
Col_SD <- function(col1,col2){
  mean <- (col1+col2)/2
  SD <- sqrt(((mean - col1)^2 + (mean - col2)^2)/2)
  return(SD)
}

### Calculate the Number of Unique Peptides Quantified ###
Num_Peptides <- function(ProtList, PepList1, PepList2){
  for(i in 1:nrow(ProtList)){
    peptides <- merge(PepList1[[i]],PepList2[[i]],all=T)
    ProtList$Unique_Peptides[i] <- length(unique(peptides$Net.peptide))
  }
  return(ProtList)
}

### Make assigned peptide list #####
Peptide_List <- function(ProtList,peptides){
  
  pullProt <- function(x){
    df <-  peptides[[x]]
    df$Assigned_Protein <- x
    df <- df[order(df$Net.peptide,df$Replicate),]
    df
  }
  AssignedList <- lapply(ProtList$Protein,pullProt)
  
  count <- sum(sapply(AssignedList,nrow))
  
  merge.all <- function(x,y){
    merge(x,y,all=T)
  }
  out <- Reduce(merge.all, AssignedList)
  out <- out[order(out$Assigned_Protein,out$Net.peptide,out$Replicate),]
  out
}

### Update Protein Groups ###
Prot_Group <- function(ProtList){
  ProtList <- ProtList[order(ProtList$Protein_Group,ProtList$Probability),]
  k <- 1
  for (i in 1:nrow(ProtList)){
    ProtList$Protein_Group[i] <- k
    if( i != nrow(ProtList) & (ProtList$Protein_Group[i] != ProtList$Protein_Group[i+1])){
      k <- k+1
    }
  }
  return(ProtList)
}