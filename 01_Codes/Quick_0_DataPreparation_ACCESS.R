# Info --------------------------------------------------------------------
#
# Author: Luca Montana
# Updated by: Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch
# Location: Institut Maurice Lamontagne
# Date: 2021-2024
#
# Overview: Prepare sequences to update haplotype libraries - Multiple Sequences Alignment (MSA)
# Multiple sequences alignment
#

# 0. Housekeeping ---------------------------------------------------------

# Libraries
# if(!require(tidyverse)){install.packages("tidyverse")}
# library(readxl)
library(dplyr)

# if(!require(data.table)){install.packages("data.table")}
library(data.table)  # rleid function

# if(!require(stringr)){install.packages("stringr")}
library(stringr)  # str_count function

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings", force = TRUE)  # install Biostrings the first time you run this script
# BiocManager::install("msa", force = TRUE)  # install msa the first time you run this script
library(Biostrings)
library(msa)

# Functions
"%nin%" <- Negate("%in%")

# 1. Data -----------------------------------------------------------------

# First give a name to this analysis (for path)

projet <- "BelugaMonitorage_Part1_20240503"


## 1.1. Upload databases --------------------------------------------------

library(RODBC)
library(remotes)

# Package maison pour faciliter l'utilisation de R pour accéder aux données
#remotes::install_github("GenomicsMLI-DFO/BD-LabGeno-IML_gestion")
library(BDLG.gestion)

test_DB()
list_DB()

data <- BDLG.gestion::load_DB("Metadata_Sequencage")

data |>  group_by(Nom_commun, Loci_sequencage) |> summarise(N = n())

dt <- data |> dplyr::filter(Nom_commun == "BELUGA", Loci_sequencage == "Dloop")

# How many without haplotype?
dt |> dplyr::filter(is.na(Haplotype)) |> group_by(Annee_sequencage) |> summarise(N = n())

dt <- dt |> dplyr::filter(Annee_sequencage == 2024)

## 1.2. Format input database for MSA -------------------------------------


### 1.2.4. Format Sequence_consensus --------------------------------------

dt$Sequence_consensus <- toupper(dt$Sequence_consensus)  # Nucleotides in capital letters

dt |> dplyr::filter(stringr::str_detect(Sequence_consensus, "-"))
#dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)  # No breaks within sequences, remove '-' on the edges

dt |> dplyr::filter(stringr::str_detect(Sequence_consensus, ":"))
#dt$Sequence_consensus <- gsub(":", "", dt$Sequence_consensus)  # No breaks within sequences, remove ':' in one sequence S_22_05102

#### 1.2.4.1. Filter by sequence length: remove short sequences -----------

dt$N_nucl <- as.integer(nchar(dt$Sequence_consensus))
dt %>% pull(N_nucl) %>% table(useNA = "ifany")
dt <- dt[!dt$N_nucl < 200,]  # remove any sequence that is very short (there should be none)
dt <- dt[!is.na(dt$Sequence_consensus),]

# Verify if and which ambiguities are included in sequences at present
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))


#### 1.2.4.2. Identify duplicated sequences -------------------------------
# Necessary step to recognize different which sequences corresponds to which extraction/re-sequencing event in a fasta file

duplicated(dt$Numero_unique_specimen) |> table()

dloop <- dt %>%  # identify duplicated sequences by adding -2, -3, -4 after the ID of the specimen
  arrange(Numero_unique_specimen, Numero_unique_extrait) %>%
  group_by(Numero_unique_specimen) %>%  # group by specimen ID
  mutate(Duplicated = row_number(Numero_unique_specimen)) %>%  # create new Duplicated column specifying which specimen is duplicated using numbers
  # previously using rleid (from data.table), but found out that some duplicated Numero_unique_specime were issued from duplicated Numero_unique_extrait
  # rleid(Numero_unique_extrait)
  mutate(Numero_unique_specimen = paste(Numero_unique_specimen, Duplicated, sep = "-"))  # paste specimen ID with duplication number (add identified to Numero_unique_specimen)
dloop$Numero_unique_specimen <- gsub("-1", "", dloop$Numero_unique_specimen)  # unnecessary to specify which ones are unique or the first specimen of a series of duplicates

nrow(dloop)

# 2. Multiple Sequence Alignment ------------------------------------------

## 2.1. Prepare DNAStringSet object ---------------------------------------

seq <- dloop$Sequence_consensus  # extract sequences
dna <- DNAStringSet(seq)  # create DNAStringSet object
names(dna) <- dloop$Numero_unique_specimen  # name each sequence in the DNAStringSet object
 writeXStringSet(dna, paste0("00_Data/01_fasta/",projet,"_complete_seq_rev_comp.fasta"))
# writeXStringSet(dna, "00_Data/01_fasta/Beluga_Narwhal_complete_seq_rev_comp.fasta")


## 2.2. Reverse complement ------------------------------------------------
# Some sequences can be in reverse complement: check out visually in MEGA if any of the new sequences is in revcomp
# Procedure: open MEGA-X, upload fasta file created above (commented out) - Align, Alignment - Align by MUSCLE, alignment options -
# MegaX (remember to add consensus sequence 615 bp at the top)
# Settings: Gap open = -10000.00; Gap Extend = -400.00; Max Iterations = 30; Cluster Method = UPGMA
# Found so far: S_20_00647 (not the duplicate S_20_00647-2), 01198, 01618, 01638, 02908, 03180, 03202

# dna$S_24_04431 <- reverseComplement( dna$S_24_04431)

# writeXStringSet(dna, "00_Data/01_fasta/Beluga_complete_seq_n3441.fasta")  # remember to change sample size if new sequences are included
# writeXStringSet(dna, "00_Data/01_fasta/Beluga_Narwhal_complete_seq_n3459.fasta")  # remember to change sample size if new sequences are included


## 2.3. MSA ---------------------------------------------------------------

dna.algn <- msa(dna, method = "Muscle", gapOpening = 10000, gapExten = 400, maxiters = 30, type = "dna",
                order = "input", verbose = T)
print(dna.algn, show = "complete")
alignment <- DNAStringSet(dna.algn)  # to save the alignment
writeXStringSet(alignment, paste0("00_Data/01_fasta/",projet, ".fasta"))  # write sample size at the end of the fasta file name
# writeXStringSet(alignment, "00_Data/01_fasta/Beluga_Narwhal_alignment_complete_n3460.fasta")  # write sample size at the end of the fasta file name
# dna.algn <- readDNAStringSet("00_Data/01_fasta/Beluga_alignment_complete_n3441.fasta")  # upload complete alignment
# dna.algn <- readDNAStringSet("00_Data/01_fasta/Beluga_Narwhal_alignment_complete_n3460.fasta")  # upload complete alignment



# 3. Standardize sequences length: short and long sequences  --------------
# Use bp strings at start and end of consensus sequences (615bp and 234bp) to define the cutting loci


## 3.1. Cut sequences - 234 bp --------------------------------------------

### 3.1.1. Define F and R 'primers' ---------------------------------------
# Lillie et al. 1996: in the 235 basepairs routinely sequenced (position 134 to 384) we have found 18 variable sites expressing 37 haplotypes in 450 individuals
# De March & Postma 2003: the mtDNA locus we used consists of 234 nucleotides that are found at the beginning of the d-loop region of the molecule (Brown Gladden et al. 1997)
# Brown-Gladden et al. 1997: first SNP = 129, last SNP = 355
# Primer 5'-3' used in Brown Gladden et al. and De March & Postma: Bel5 <- "ACATTTTACTGTGACTATTG"  # at the beginning of 5'-end of CR (nt = 71)

cr917 <- readDNAStringSet(filepath = "00_Data/01_fasta/Sequence_Dloop_complete.fasta")  # NCBI ID: U18117.1
F234 <- toString(subseq(cr917, start = 126, width = 21))  # at the beginning of 5'-end of CR (nt = 126)
R234 <- toString(subseq(cr917, end = 359, width = 21))


### 3.1.2. Define cutting position of alignment ----------------------------

res.F234 <- vmatchPattern(DNAString(F234), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F234.int <- res.F234@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.F234 <- as.numeric(as.character(cut.F234.int[[1]])) - nchar(F234) + 1  # position 175 - if multiple position it means at least one sequence was misaligned
cut.F234

res.R234 <- vmatchPattern(DNAString(R234), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R234.int <- res.R234@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.R234 <- as.numeric(as.character(cut.R234.int[[1]]))  # position 408 - if multiple position it means at least one sequence was misaligned
cut.R234

### 3.1.3. Cut sequences - save fasta --------------------------------------

Dloop234 <- subseq(DNAStringSet(dna.algn), start = cut.F234, end = cut.R234)
print(Dloop234, show = "complete")
table(Dloop234@ranges@width)
writeXStringSet(Dloop234, paste0("00_Data/01_fasta/", projet, "_234bp_n", length(Dloop234), ".fasta"))  # save fasta
#Dloop234 <- readDNAStringSet("fasta/Beluga_234bp_n3441.fasta")  # upload 234bp alignment


### 3.1.4. Save dataset -----------------------------------------------------

dna234 <- data.frame(ID = names(Dloop234),
                     Sequence = Dloop234) %>%
  left_join(dloop[,c("Numero_unique_specimen","Numero_unique_extrait",
                     "Numero_unique_groupe", "Region_echantillonnage", "Lieu_echantillonnage", "Annee_echantillonnage", "Mois_echantillonnage", "Jour_echantillonnage",
                     "Amorce_sequencage_F", "Amorce_sequencage_R", "Annee_sequencage", "Mois_sequencage"  )], by = c("ID"="Numero_unique_specimen"))
write.table(dna234, file = paste0("00_Data/02_dloop_clean/", projet, "_234bp_n", length(Dloop234),"_metadata.txt"), row.names = F)

write.table(dna234, file = paste0("00_Data/02_dloop_clean/", projet, "_metadata.txt"), row.names = F)

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
writeXStringSet(Dloop570,  paste0("00_Data/01_fasta/", projet, "_570bp_n", length(Dloop570), ".fasta"))  # save fasta
#Dloop615 <- readDNAStringSet("fasta/Beluga_615bp_n3441.fasta")  # upload 615bp alignment


### 3.1.4. Save dataset -----------------------------------------------------

dna570 <- data.frame(ID = names(Dloop570),
                     Sequence = Dloop570) %>%
  left_join(dloop[,c("Numero_unique_specimen","Numero_unique_extrait",
                     "Numero_unique_groupe", "Region_echantillonnage", "Lieu_echantillonnage", "Annee_echantillonnage", "Mois_echantillonnage", "Jour_echantillonnage",
                     "Amorce_sequencage_F", "Amorce_sequencage_R", "Annee_sequencage", "Mois_sequencage"  )], by = c("ID"="Numero_unique_specimen"))
write.table(dna570, file = paste0("00_Data/02_dloop_clean/", projet, "_570bp_n", length(Dloop570),"_metadata.txt"), row.names = F)

## 3.3. Cut sequences - 615 bp --------------------------------------------

### 3.3.1. Define F and R 'primers' ---------------------------------------
# Use 5'- and 3'- ends of consensus sequence

cr615 <- readDNAStringSet(filepath = "00_Data/01_fasta/Sequence_Dloop_Ref_615pb_consensus.fasta")
F615 <- toString(subseq(cr615, start = 1, width = 21))  # at the beginning of 5'-end of CR (nt = 38)
R615 <- toString(subseq(cr615, end = 615, width = 21))


### 3.3.2. Define cutting position of alignment ----------------------------

res.F615 <- vmatchPattern(DNAString(F615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F615.int <- res.F615@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.F615 <- as.integer(as.character(cut.F615.int[[1]])) - nchar(F615) + 1  # position 87 - if multiple position it means at least one sequence was misaligned
cut.F615

res.R615 <- vmatchPattern(DNAString(R615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R615.int <- res.R615@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.R615 <- as.numeric(as.character(cut.R615.int[[1]]))  # position 701 - if multiple position it means at least one sequence was misaligned
cut.R615

### 3.3.3. Cut sequences - save fasta --------------------------------------

Dloop615 <- subseq(DNAStringSet(dna.algn), start = cut.F615, end = cut.R615)
print(Dloop615, show = "complete")
table(Dloop615@ranges@width)
Dloop615
writeXStringSet(Dloop615,  paste0("00_Data/01_fasta/", projet, "_615bp_n", length(Dloop615), ".fasta"))  # save fasta
#Dloop615 <- readDNAStringSet("fasta/Beluga_615bp_n3441.fasta")  # upload 615bp alignment


### 3.3.4. Save dataset -----------------------------------------------------

dna615 <- data.frame(ID = names(Dloop615),
                     Sequence = Dloop615) %>%
  left_join(dloop[,c("Numero_unique_specimen","Numero_unique_extrait",
                     "Numero_unique_groupe", "Region_echantillonnage", "Lieu_echantillonnage", "Annee_echantillonnage", "Mois_echantillonnage", "Jour_echantillonnage",
                     "Amorce_sequencage_F", "Amorce_sequencage_R", "Annee_sequencage", "Mois_sequencage"  )], by = c("ID"="Numero_unique_specimen"))
write.table(dna615, file = paste0("00_Data/02_dloop_clean/", projet, "_615bp_n", length(Dloop615),"_metadata.txt"), row.names = F)


# 4. Prepare datasets to run Script1 --------------------------------------

nt <- c("A","T","C","G")  # define DNA bases
ambiguous <- c("N","R","Y","K","M","S","W","B","D","H","V") # define ambiguities


## 4.1. Short sequences ---------------------------------------------------

### 4.1.1. Info for dataset -----------------------------------------------
# Include number of nt, number of ACTG, number ambiguities and number of missing nt in dataset

seq234 <- dna234$Sequence
exp_seq_len234 <- 234

info234 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info234) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq234)){
  if(is.na(seq234[i])){
    info234[i,] <- c("NA","NA","NA","NA")
  }else{
    seq_len <- sum(str_count(seq234[i], c(nt, ambiguous)))
  }
  info234[i,] <- c(seq_len, sum(str_count(seq234[i], nt)), sum(str_count(seq234[i], ambiguous)), exp_seq_len234-seq_len)
}
dna234 <- cbind(dna234, info234)

### 4.1.2. Filter dataset ---------------------------------------------------

# Remove sequences with ambiguities or with less than 234 nt
rem <- dna234[dna234$N.ATCG < 234, "ID"]
dna234_red <- dna234[!(dna234$ID %in% rem),]

# Remove duplicated sequences
dup <- dna234_red[base::grepl("-", dna234_red$ID), "ID"]
dna234_red <- dna234_red[!(dna234_red$ID %in% dup), ]


### 4.1.3. Save clean fasta -------------------------------------------------

seq234_red <- dna234_red$Sequence
seq234_red
s234.red <- DNAStringSet(seq234_red)
names(s234.red) <- dna234_red$ID

length(seq234_red)

writeXStringSet(s234.red, paste0("00_Data/01_fasta/", projet, "_234bp_onlyATGC_n",length(seq234_red), ".fasta"))


## 4.2. Long sequences ----------------------------------------------------

### 4.2.1. Info for dataset -----------------------------------------------
# Include number of nt, number of ACTG, number ambiguities and number of missing nt in dataset

seq570 <- dna570$Sequence
exp_seq_len570 <- 570

info570 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info570) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq570)){
    if(is.na(seq570[i])){
        info570[i,] <- c("NA","NA","NA","NA")
    }else{
        seq_len <- sum(str_count(seq570[i],c(nt, ambiguous)))
    }
    info570[i,] <- c(seq_len, sum(str_count(seq570[i], nt)), sum(str_count(seq570[i], ambiguous)), exp_seq_len570-seq_len)
}
dna570 <- cbind(dna570, info570)


### 4.2.2. Filter dataset ---------------------------------------------------

# Remove sequences with ambiguities or with less than 234 nt
rem <- dna570[dna570$N.ATCG < 570, "ID"]
dna570_red <- dna570[!(dna570$ID %in% rem),]

# Remove duplicated sequences
dup <- dna570_red[base::grepl("-", dna570_red$ID), "ID"]
dna570_red <- dna570_red[!(dna570_red$ID %in% dup), ]


dna570[(dna570$ID %in% rem),] |> dplyr::select(ID, N.nucl, N.ATCG, N.ambig, N.manquants)

dna570 |> dplyr::filter(ID == "S_24_04858")

### 4.1.3. Save clean fasta -------------------------------------------------

seq570_red <- dna570_red$Sequence
s570.red <- DNAStringSet(seq570_red)
names(s570.red) <- dna570_red$ID

length(s570.red)
s570.red

writeXStringSet(s570.red,  paste0("00_Data/01_fasta/", projet, "_570bp_onlyATGC_n",length(seq570_red), ".fasta"))

