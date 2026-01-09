# Info --------------------------------------------------------------------
#
# Author: Luca Montana
# Updated by: Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch
# Location: Institut Maurice Lamontagne
# Date: 2021-2026
#
# Overview: Prepare sequences to update haplotype libraries - Multiple Sequences Alignment (MSA)
# Multiple sequences alignment
#

# 0. Housekeeping ---------------------------------------------------------

# Libraries
if(!require(dplyr)){install.packages("dplyr")}
if(!require(data.table)){install.packages("data.table")}
if(!require(stringr)){install.packages("stringr")}
if(!require(readr)){install.packages("readr")}

library(dplyr)
library(data.table)  # rleid function
library(stringr)  # str_count function
library(readr)

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings", force = TRUE)  # install Biostrings the first time you run this script
# BiocManager::install("msa", force = TRUE)  # install msa the first time you run this script
library(Biostrings)
library(msa)

# Package maison pour faciliter l'utilisation de R pour accéder aux données
#remotes::install_github("GenomicsMLI-DFO/BD-LabGeno-IML_gestion")
library(BDLG.gestion)

# Functions
"%nin%" <- Negate("%in%")

# 1. Data -----------------------------------------------------------------

# First give a name to this analysis (for path)

projet <- "Monitorage_20260108"

if (!dir.exists(file.path("00_Data/02_dloop_clean/", projet)))
  dir.create(file.path("00_Data/02_dloop_clean/", projet))

## 1.1. Upload databases from ACCESS repository ---------------------------


test_DB()
list_DB()

dt <- BDLG.gestion::load_DB("R_Metadata_Sequencage_Dloop")

# How many without haplotype?
dt |> dplyr::filter(is.na(Haplotype)) |> group_by(Annee_sequencage) |> summarise(N = n())

dt <- dt |> dplyr::filter((Annee_sequencage == 2025 & Mois_sequencage == 12) | Annee_sequencage == 2026)

nrow(dt)

## 1.2. Format input database for MSA -------------------------------------

source("./01_Codes/R_functions.R")

verif_consensus_sequences(dt, seq_col = "Sequence_consensus", min_length = 570 )
# Si NULL/message, tout est OK
# Sinon, audit contient le tableau des anomalies avec leurs positions pour corriger efficacement

# Pour éliminer les séquences NA
dloop <- dt |> dplyr::filter(!is.na(Sequence_consensus))

verif_consensus_sequences(dloop, seq_col = "Sequence_consensus", min_length = 570 )


# 2. Multiple Sequence Alignment ------------------------------------------

## 2.1. Prepare DNAStringSet object ---------------------------------------

seq <- dloop$Sequence_consensus  # extract sequences
dna <- DNAStringSet(seq)  # create DNAStringSet object
names(dna) <- dloop$ID # name each sequence in the DNAStringSet object

#writeXStringSet(dna, paste0("00_Data/01_fasta/",projet,"_complete_seq_rev_comp.fasta"))

## 2.3. MSA ---------------------------------------------------------------

dna.algn <- msa(dna, method = "Muscle", gapOpening = 10000, gapExten = 400, maxiters = 30, type = "dna",
                order = "input", verbose = T)
print(dna.algn, show = "complete")
alignment <- DNAStringSet(dna.algn)  # to save the alignment

#writeXStringSet(alignment, paste0("00_Data/01_fasta/",projet, ".fasta"))  # write sample size at the end of the fasta file name

# 3. Standardize sequences length: short and long sequences  --------------
## 3.2. Cut sequences - 570 bp --------------------------------------------

### 3.2.1. Define F and R 'primers' ---------------------------------------
# Use 5'- and 3'- ends of consensus sequence

cr570 <- readDNAStringSet(filepath = "00_Data/01_fasta/Sequence_Dloop_Ref_570pb_minimale.fasta")
F570 <- toString(subseq(cr570, start = 1, width = 21))  # at the beginning of 5'-end of CR (nt = 38)
R570 <- toString(subseq(cr570, end = 570, width = 21))


### 3.2.2. Define cutting position of alignment ----------------------------

res.F570 <- vmatchPattern(DNAString(F570), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F570.int <- res.F570@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.F570 <- as.integer(as.character(cut.F570.int[[1]])) - nchar(F570) + 1  # position 87 - if multiple position it means at least one sequence was misaligned
cut.F570

res.R570 <- vmatchPattern(DNAString(R570), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R570.int <- res.R570@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.R570 <- as.numeric(as.character(cut.R570.int[[1]]))  # position 701 - if multiple position it means at least one sequence was misaligned
cut.R570

### 3.2.3. Cut sequences - save fasta --------------------------------------

Dloop570 <- subseq(DNAStringSet(dna.algn), start = cut.F570, end = cut.R570)
print(Dloop570, show = "complete")
table(Dloop570@ranges@width)
Dloop570

writeXStringSet(Dloop570,  paste0("00_Data/02_dloop_clean/", projet, "/", projet, "_570bp.fasta"))  # save fasta


### 3.1.4. Save dataset -----------------------------------------------------

dna570 <- data.frame(ID = names(Dloop570),
                     Sequence = Dloop570) %>%
  full_join(dt |> mutate(ID = as.character(ID)))
readr::write_csv(dna570, file = paste0("00_Data/02_dloop_clean/", projet, "/", projet, "_metadata.csv"))

verif_consensus_sequences(dna570, seq_col = "Sequence", min_length = 570 )


