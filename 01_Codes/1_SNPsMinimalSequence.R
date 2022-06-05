# Info --------------------------------------------------------------------
# 
# Author: Benjamin Hornoy
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-12-17
# 
# Overview: Detection of variation and minimal sequences
# WARNING: ONLY complete sequences (615 or 234 bp) and without ambiguities or gaps should be used
# See file: "Pipeline bélugas - Benjamin Hornoy - 13082021.pdf"
# 
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
# if(!require(adegenet)){install.packages("adegenet")}
library(adegenet)

# Functions




# 1. Data -----------------------------------------------------------------
# Use ADEGENET to import fasta

myDNA615 <- fasta2DNAbin("00_Data/01_fasta/Beluga_615bp_onlyATGC_n3329.fasta")
seq_len615 <- 615  # length of expected sequence
myDNA234 <- fasta2DNAbin("00_Data/01_fasta/Beluga_234bp_onlyATGC_n3423.fasta")
seq_len234 <- 234  # length of expected sequence




# 2. Identify polymorphic sites -------------------------------------------
# Read sequences in ADEGENET and specify expected length of sequence
# Conserve only polymorphic sites

## 2.1. Short sequences ---------------------------------------------------
# Brown Gladden et al. 1997 found 19 polymorphic sites (minimal seq = 226 nt)

obj234 <- DNAbin2genind(myDNA234, polyThres=0)
obj234  # for information
snpPos234 <- locNames(obj234)


## 2.2. Long sequences ----------------------------------------------------

obj615 <- DNAbin2genind(myDNA615, polyThres=0)  # polyThres defines the minimum frequency of an allele to considerate locus as polymorphic
obj615  # for information
snpPos615 <- locNames(obj615)




# 3. Summary of results ---------------------------------------------------

## 3.1. Short sequences ---------------------------------------------------

snps234 <- data.frame()
snps234 <- cbind(length(snpPos234), seq_len234, length(snpPos234)/seq_len234*100, paste(snpPos234[1],"-",snpPos234[length(snpPos234)]),
                 as.numeric(snpPos234[length(snpPos234)])-as.numeric(snpPos234[1])+1, paste(snpPos234, collapse=","), length(obj234@ploidy))
colnames(snps234) <- c("Nb.de.SNPs", "Longueur.de.la.sequence", "%.de.snps", "Bornes.sequence.minimale", "Longueur.sequence.minimale",
                       "Position.des.SNPS", "Taille.echantillons")


## 3.2. Long sequences ----------------------------------------------------

snps615 <- data.frame()
snps615 <- cbind(length(snpPos615), seq_len615, length(snpPos615)/seq_len615*100, paste(snpPos615[1],"-",snpPos615[length(snpPos615)]),
                 as.numeric(snpPos615[length(snpPos615)])-as.numeric(snpPos615[1])+1, paste(snpPos615, collapse=","), length(obj615@ploidy))
colnames(snps615) <- c("Nb.de.SNPs", "Longueur.de.la.sequence", "%.de.snps", "Bornes.sequence.minimale", "Longueur.sequence.minimale",
                       "Position.des.SNPS", "Taille.echantillons")



# 4. Save dataset ---------------------------------------------------------

snps <- rbind(snps615, snps234)  # merge tables
write.csv(snps, file = "02_Results/01_poly_seq_min/polymorphismes_et_seq_minimale.csv", row.names = F)

