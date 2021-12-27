# Prepare sequence dataset for haplotype assignment:
# Detection of variation and minimal sequences
# 
# Benjamin Hornoy
# 
#

getwd()
rm(list = ls())

# Libraries and function --------------------------------------------------
if(!require(adegenet)){install.packages("adegenet")}
library(adegenet)


# Data --------------------------------------------------------------------
# Use ADEGENET to import fasta
myDNA <- fasta2DNAbin("Beluga_615bp_n3284.fasta")
seq_length <- 615 #longueur de la séquence attendue


# Identify polymorphic sites ----------------------------------------------
# Read sequences in ADEGENET and specify expected length of sequence
# Conserve only polymorphic sites
obj <- DNAbin2genind(myDNA, polyThres=0)
obj #pour information
snpPos <- locNames(obj)

# Assign summary of results in a table
snps <- data.frame()
snps <- cbind(length(snpPos), seq_length, length(snpPos)/seq_length*100, paste(snpPos[1],"-",snpPos[length(snpPos)]), as.numeric(snpPos[length(snpPos)])-as.numeric(snpPos[1])+1, paste(snpPos, collapse=","))
colnames(snps) <- c("Nb.de.SNPs", "Longueur.de.la.sequence", "%.de.snps", "Bornes.sequence.minimale", "Longueur.sequence.minimale", "Position.des.SNPS")
write.csv(snps, file="polymorphismes_et_seq_minimale_n3284.csv", row.names = F)


