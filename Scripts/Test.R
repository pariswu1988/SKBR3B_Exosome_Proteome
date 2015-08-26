#! c:\"Program Files"\R\R-3.1.1\bin\x64\Rscript --vanilla

Match_ID <- function(x,field) {
 b <- Proteins[x == Proteins$Protein_ID,field]
 return(b)
}

test <- sapply(SI2$Assigned_Protein, Match_ID,"Description")
