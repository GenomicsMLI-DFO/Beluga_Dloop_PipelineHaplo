# Info --------------------------------------------------------------------
#
# Authors: Benjamin Hornoy, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch
# Location: Institut Maurice Lamontagne
# Date: 2021-12-17
#
# Overview: Assign haplotypes to specimens
#
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
#getwd()

# Libraries
#library(readxl)
library(dplyr)

# Functions
"%nin%" <- Negate("%in%")


# 1. Data -----------------------------------------------------------------

data234 <- read.table("00_Data/02_dloop_clean/Sequences_Dloop234_n14_Archambault.txt", header = T)
str(data234)
data615 <- read.table("00_Data/02_dloop_clean/Sequences_Dloop615_n14_Archambault.txt", header = T)
str(data615)

# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in 2a_HaploLibrary_234.R and 2b_HaploLibrary_615.R

data615 |> dplyr::select(Numero_un)

s234.red <- Biostrings::readDNAStringSet("00_Data/01_fasta/Beluga_234bp_onlyATGC_n14_Archambault.fasta")

s615.red <- Biostrings::readDNAStringSet("00_Data/01_fasta/Beluga_615bp_onlyATGC_n14_Archambault.fasta")

s234.df <- data.frame(Numero_unique_specimen = names(s234.red),
                   seq = as.character(s234.red)) |> dplyr::left_join(lib234) |> dplyr::select(-seq)

s615.df <- data.frame(Numero_unique_specimen = names(s615.red),
                      seq = as.character(s615.red)) |> dplyr::left_join(lib615) |> dplyr::select(-seq)

s234.df |> dplyr::full_join(s615.df)

# 2. Assign haplotype to each individual ----------------------------------

## 2.1. Upload haplotype libraries ------------------------------------------

lib234 <- read.csv('02_Results/00_libraries/librairie_54_haplotypes234.csv')  # most recent haplotype library
colnames(lib234) <- c("hapl_234","seq")  # if it's not already the case
lib615 <- read.csv('02_Results/00_libraries/librairie_143_haplotypes615.csv')  # most recent haplotype library
colnames(lib615) <- c("hapl_615","seq")  # if it's not already the case


