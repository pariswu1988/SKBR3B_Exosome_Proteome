carta <- read.delim("Data/exocarta.txt")
prot <- read.csv("Temp/Proteins.csv", stringsAsFactors = F)
prot <- prot[prot$Unique_Peptides >= 3,]
prot$prot.names <- prot$Protein
prot$prot.names <- gsub("^.*\\|","", prot$prot.names)
prot$prot.names <- gsub("_.*$","",prot$prot.names)

m <- pRolocmarkers("hsap")
m <- data.frame("prot.names" = names(m), "localization" = m)
m$prot.names <- gsub("_.*$","", m$prot.names)

m2 <- m[m$prot.names %in% prot$prot.names,]
m2

prot2 <- merge(prot,m2)

carta <- unique(as.character(carta[carta$SPECIES == "Homo sapiens", "GENE.SYMBOL"]))

protMarkersTop <- prot2[!(prot2$Gene_Name %in% carta),]
protMarkersTop

others <- prot[!(prot$Gene_Name %in% carta) & prot$Probability == 1,]
others <- others[order(others$Unique_Peptides, decreasing = T),]

picked <- c("HK2", "TMX3", "MST4", "SEC62", "TECR")
localizations <- c("Mitochondrion","ER", "Cytoplasm", "ER","ER")

others <- others[others$Gene_Name %in% picked,]
others$localization <- localizations

nonExoMarkers <- merge(others, protMarkersTop, all = T)
nonExoMarkers$Localization <- rep("Non-Exosome", nrow(nonExoMarkers))

write.csv(nonExoMarkers[,c("Protein","Localization")], "Temp/nonExosomeMarkers.csv",
          row.names=F)
