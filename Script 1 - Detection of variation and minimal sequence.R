# Prepare sequence dataset for haplotype assignment:
# Detection of variation and minimal sequences
# WARNING: ONLY complete sequences (615 or 234 bp) and without ambiguities or gaps should be used
# See Pipeline bélugas - Benjamin Hornoy - 13082021.pdf
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
myDNA615 <- fasta2DNAbin("Beluga_615bp_onlyATGC_n3107.fasta")
seq_len615 <- 615  # length of expected sequence
myDNA234 <- fasta2DNAbin("Beluga_234bp_onlyATGC_n3176.fasta")
seq_len234 <- 234  # length of expected sequence


# Identify polymorphic sites ----------------------------------------------
# Read sequences in ADEGENET and specify expected length of sequence
# Conserve only polymorphic sites
# 615 bp
obj615 <- DNAbin2genind(myDNA615, polyThres=0)
obj615  # for information
snpPos615 <- locNames(obj615)

# 234 bp
# Brown Gladden et al. 1997 found 19 polymorphic sites (minimal seq = 226 nt)
obj234 <- DNAbin2genind(myDNA234, polyThres=0)
obj234  # for information
snpPos234 <- locNames(obj234)


# Assign summary of results in a table ------------------------------------
# 615 bp
snps615 <- data.frame()
snps615 <- cbind(length(snpPos615), seq_len615, length(snpPos615)/seq_len615*100, paste(snpPos615[1],"-",snpPos615[length(snpPos615)]),
                 as.numeric(snpPos615[length(snpPos615)])-as.numeric(snpPos615[1])+1, paste(snpPos615, collapse=","), length(obj615@ploidy))
colnames(snps615) <- c("Nb.de.SNPs", "Longueur.de.la.sequence", "%.de.snps", "Bornes.sequence.minimale", "Longueur.sequence.minimale",
                       "Position.des.SNPS", "Taille.echantillons")

# 234 bp
snps234 <- data.frame()
snps234 <- cbind(length(snpPos234), seq_len234, length(snpPos234)/seq_len234*100, paste(snpPos234[1],"-",snpPos234[length(snpPos234)]),
                 as.numeric(snpPos234[length(snpPos234)])-as.numeric(snpPos234[1])+1, paste(snpPos234, collapse=","), length(obj234@ploidy))
colnames(snps234) <- c("Nb.de.SNPs", "Longueur.de.la.sequence", "%.de.snps", "Bornes.sequence.minimale", "Longueur.sequence.minimale",
                       "Position.des.SNPS", "Taille.echantillons")


# Merge and save dataset --------------------------------------------------
snps <- rbind(snps615, snps234)
write.csv(snps, file = "polymorphismes_et_seq_minimale.csv", row.names = F)

