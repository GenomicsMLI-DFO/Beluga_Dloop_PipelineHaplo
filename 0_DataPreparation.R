# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-12-17
# 
# Overview: Prepare sequences to update haplotype libraries - Multiple Sequences Alignment (MSA)
# Multiple sequences alignment
# 
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
# if(!require(tidyverse)){install.packages("tidyverse")}
library(readxl)
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

## 1.1. Upload databases --------------------------------------------------

# Originally in ACCESS folder on Drive. Specify the path to the directory where the file is stored
d <- read_excel("../../MOBELS/DB/ACCESS/20220524_MOBELS.xlsx", sheet = "D-Loop", na = "NA")  # remember to specify right path to beluga ACCESS dataset
s <- read_excel("../../MOBELS/DB/ACCESS/20220524_MOBELS.xlsx", sheet = "Specimens", na = "NA")  # remember to specify right path to beluga ACCESS dataset
g <- read_excel("../../MOBELS/DB/ACCESS/20220524_MOBELS.xlsx", sheet = "Groupe", na = "NA")  # remember to specify right path to beluga ACCESS dataset 

# d <- read.csv("Dloop_MOBELS.csv")
# s <- read_excel("../ACCESS/20220328_MOBELS.xlsx", sheet = "Specimens", na = "NA")



## 1.2. Format input database for MSA -------------------------------------

### 1.2.1. Dloop ----------------------------------------------------------

str(d)  # 3993 rows
colnames(d)[2] <- "Numero_unique_extrait"

# Subset dataset: remove 'useless' columns
# d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
# d <- transform(d, Qualite_sequence = as.integer(d$Qualite_sequence),
#                N_nucl = as.integer(d$N_nucl))
d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Sequence_consensus))

# Remove specimens without consensus sequence
d <- d[!is.na(d$Sequence_consensus),]  # removes 469 rows


### 1.2.2. Specimens ------------------------------------------------------

str(s)

# Subset datase: remove 'useless' columns
s <- s[,c("Numero_unique_specimen","Nom_commun")]


### 1.2.3. Merge 'd' and 's' datasets -------------------------------------

dt <- merge(d, s, by = "Numero_unique_specimen")


#### 1.2.3.1. Species: remove narwhals and putative hybrids ---------------

table(dt$Nom_commun, useNA = 'ifany')  # 3510 beluga; 2 hybrids; 12 narwhals
dt <- dt[dt$Nom_commun %in% "Beluga", colnames(dt) %nin% "Nom_commun"]  # removes 14 specimens and Nom_commun column


### 1.2.4. Format Sequence_consensus --------------------------------------

dt$Sequence_consensus <- toupper(dt$Sequence_consensus)  # Nucleotides in capital letters
dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)  # No breaks within sequences, remove '-' on the edges


#### 1.2.4.1. Filter by sequence length: remove short sequences -----------

dt$N_nucl <- as.integer(nchar(dt$Sequence_consensus))
table(dt$N_nucl)
dt <- dt[!dt$N_nucl < 200,]  # remove any sequence that is very short (there should be none)

# Verify if and which ambiguities are included in sequences at present
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))


#### 1.2.4.2. Identify duplicated sequences -------------------------------
# Necessary step to recognize different which sequences corresponds to which extraction/re-sequencing event in a fasta file

dloop <- dt %>%  # identify duplicated sequences by adding -2, -3, -4 after the ID of the specimen
  group_by(Numero_unique_specimen) %>%  # group by specimen ID
  mutate(Duplicated = rleid(Numero_unique_extrait)) %>%  # create new Duplicated column specifying which specimen is duplicated using numbers
  mutate(Numero_unique_specimen = paste(Numero_unique_specimen, Duplicated, sep = "-"))  # paste specimen ID with duplication number (add identified to Numero_unique_specimen)
dloop$Numero_unique_specimen <- gsub("-1", "", dloop$Numero_unique_specimen)  # unnecessary to specify which ones are unique or the first specimen of a series of duplicates




# 2. Multiple Sequence Alignment ------------------------------------------

## 2.1. Prepare DNAStringSet object ---------------------------------------

seq <- dloop$Sequence_consensus  # extract sequences
dna <- DNAStringSet(seq)  # create DNAStringSet object
names(dna) <- dloop$Numero_unique_specimen  # name each sequence in the DNAStringSet object
# writeXStringSet(dna, "fasta/Beluga_complete_seq_rev_comp.fasta")


## 2.2. Reverse complement ------------------------------------------------
# Some sequences can be in reverse complement: check out visually in MEGA if any of the new sequences is in revcomp
# Procedure: open MEGA-X, upload fasta file created above (commented out) - Align, Alignment - Align by MUSCLE, alignment options -
# MegaX (remember to add consensus sequence 615 bp at the top)
# Settings: Gap open = -10000.00; Gap Extend = -400.00; Max Iterations = 30; Cluster Method = UPGMA
# Found so far: S_20_00647 (not the duplicate S_20_00647-2), 01198, 01618, 01638, 02908, 03180, 03202

dna$S_20_00647 <- reverseComplement(dna$S_20_00647)
dna$S_20_01198 <- reverseComplement(dna$S_20_01198)
dna$S_20_01618 <- reverseComplement(dna$S_20_01618)
dna$S_20_01638 <- reverseComplement(dna$S_20_01638)
dna$S_20_02908 <- reverseComplement(dna$S_20_02908)
dna$S_20_03180 <- reverseComplement(dna$S_20_03180)
dna$S_20_03202 <- reverseComplement(dna$S_20_03202)
dna$S_22_05078 <- reverseComplement(dna$S_22_05078)
dna$S_22_05196 <- reverseComplement(dna$S_22_05196)  # F surely to be sequenced
dna$S_22_05208 <- reverseComplement(dna$S_22_05208)  # F surely to be sequenced
# writeXStringSet(dna, "Beluga_complete_seq_n3510.fasta")  # remember to change sample size if new sequences are included
# weird new sequences: S_22_05057 seems fairly good but something if clearly off (maybe missing Rev? Ask Claudie about its quality);
# S_22_05080: one nt too many at pos nt 209 (either G or A if G is a mutation), also another G too many about 10--11 nt later, same as G 10-11 later again
# S_22_05081: possibly missing a C in pos 583
# S_22_05127: A G too many (nt 135 of sequence)
# S_22_05130: a CA too many (nt 481-482 of sequence);
# S_22_05133: surely missing the rev as sequence is short and end is 'crappy';
# S_22_05141: possibly a G too many at the beginning of the sequence? (64th nt);
# S_22_05160: possibly a T (or a G in case of a mutation) too many around pos 645/649 (including within seq --)
# S_22_05170: surely missing a nt aroung nt 44 of seq
# S_22_05172: same nt missing (nt 30 of seq in this case)
# S_22_05209: one A too many about pos 331 (including spaced in MEGA)





## 2.3. MSA ---------------------------------------------------------------

dna.algn <- msa(dna, method = "Muscle", gapOpening = 10000, gapExten = 400, maxiters = 30, type = "dna",
                order = "input", verbose = T)
print(dna.algn, show = "complete")
alignment <- DNAStringSet(dna.algn)  # to save the alignment
writeXStringSet(alignment, "fasta/Beluga_alignment_complete_n3314.fasta")  # write sample size at the end of the fasta file name
#dna.algn <- readDNAStringSet("Beluga_alignment_complete_n3314.fasta")  # upload complete alignment



# 3. Standardize sequences length: short and long sequences  --------------
# Use bp strings at start and end of consensus sequences (615bp and 234bp) to define the cutting loci


## 3.1. Cut sequences - 234 bp --------------------------------------------

### 3.1.1. Define F and R 'primers' ---------------------------------------
# Lillie et al. 1996: in the 235 basepairs routinely sequenced (position 134 to 384) we have found 18 variable sites expressing 37 haplotypes in 450 individuals
# De March & Postma 2003: the mtDNA locus we used consists of 234 nucleotides that are found at the beginning of the d-loop region of the molecule (Brown Gladden et al. 1997)
# Brown-Gladden et al. 1997: first SNP = 129, last SNP = 355
# Primer 5'-3' used in Brown Gladden et al. and De March & Postma: Bel5 <- "ACATTTTACTGTGACTATTG"  # at the beginning of 5'-end of CR (nt = 71)

cr917 <- readDNAStringSet(filepath = "fasta/Sequence_Dloop_complete.fasta")  # NCBI ID: U18117.1
F234 <- toString(subseq(cr917, start = 126, width = 21))  # at the beginning of 5'-end of CR (nt = 126)
R234 <- toString(subseq(cr917, end = 359, width = 21))


### 3.1.2. Define cutting position of alignment ----------------------------

res.F234 <- vmatchPattern(DNAString(F234), DNAStringSet(dna.algn), max.mismatch = 2)
cut.F234.int <- res.F234@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.F234 <- as.numeric(as.character(cut.F234.int[[1]])) - nchar(F234) + 1  # position 175

res.R234 <- vmatchPattern(DNAString(R234), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R234.int <- res.R234@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.R234 <- as.numeric(as.character(cut.R234.int[[1]]))  # position 408


### 3.1.3. Cut sequences - save fasta --------------------------------------

Dloop234 <- subseq(DNAStringSet(dna.algn), start = cut.F234, end = cut.R234)
print(Dloop234, show = "complete")
table(Dloop234@ranges@width)
writeXStringSet(Dloop234, "fasta/Beluga_234bp_n3314.fasta")  # save fasta
#Dloop234 <- readDNAStringSet("fasta/Beluga_234bp_n3314.fasta")  # upload 234bp alignment


### 3.1.4. Save dataset -----------------------------------------------------

dna234 <- data.frame(ID = names(Dloop234),
                     Sequence = Dloop234)
dna234 <- left_join(dna234, dloop[,c("Numero_unique_specimen","Numero_unique_extrait")], by = c("ID"="Numero_unique_specimen"))
write.table(dna234, file = "Sequences_Dloop234_n3314.txt", row.names = F)


## 3.2. Cut sequences - 615 bp --------------------------------------------

### 3.2.1. Define F and R 'primers' ---------------------------------------
# Use 5'- and 3'- ends of consensus sequence

cr615 <- readDNAStringSet(filepath = "fasta/Sequence_Dloop_Ref_615pb_consensus.fasta")
F615 <- toString(subseq(cr615, start = 1, width = 21))  # at the beginning of 5'-end of CR (nt = 38)
R615 <- toString(subseq(cr615, end = 615, width = 21))


### 3.1.2. Define cutting position of alignment ----------------------------

res.F615 <- vmatchPattern(DNAString(F615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F615.int <- res.F615@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.F615 <- as.integer(as.character(cut.F615.int[[1]])) - nchar(F615) + 1  # position 87

res.R615 <- vmatchPattern(DNAString(R615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R615.int <- res.R615@ends %>%
  unlist() %>%
  table() %>%
  dimnames(.)
cut.R615 <- as.numeric(as.character(cut.R615.int[[1]]))  # position 701


### 3.1.3. Cut sequences - save fasta --------------------------------------

Dloop615 <- subseq(DNAStringSet(dna.algn), start = cut.F615, end = cut.R615)
print(Dloop615, show = "complete")
table(Dloop615@ranges@width)
Dloop615
writeXStringSet(Dloop615, "fasta/Beluga_615bp_n3314.fasta")  # save fasta
#Dloop615 <- readDNAStringSet("Beluga_615bp_n3314.fasta")  # upload 615bp alignment


### 3.1.4. Save dataset -----------------------------------------------------

dna615 <- data.frame(ID = names(Dloop615),
                     Sequence = Dloop615)
dna615 <- left_join(dna615, dloop[,c("Numero_unique_specimen","Numero_unique_extrait")], by = c("ID"="Numero_unique_specimen"))
write.table(dna615, file = "Sequences_Dloop615_n3314.txt", row.names = F)




# 4. Prepare datasets to run Script1 --------------------------------------

nt <- c("A","T","C","G")  # define DNA bases
ambiguous <- c("N","R","Y","K","M","S","W","B","D","H","V") # define abmbiguities


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
s234.red <- DNAStringSet(seq234_red)
names(s234.red) <- dna234_red$ID
writeXStringSet(s234.red, "fasta/Beluga_234bp_onlyATGC_n3175.fasta")


## 4.2. Long sequences ----------------------------------------------------

### 4.2.1. Info for dataset -----------------------------------------------
# Include number of nt, number of ACTG, number ambiguities and number of missing nt in dataset

seq615 <- dna615$Sequence
exp_seq_len615 <- 615

info615 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info615) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq615)){
    if(is.na(seq615[i])){
        info615[i,] <- c("NA","NA","NA","NA")	
    }else{
        seq_len <- sum(str_count(seq615[i],c(nt, ambiguous)))
    }	
    info615[i,] <- c(seq_len, sum(str_count(seq615[i], nt)), sum(str_count(seq615[i], ambiguous)), exp_seq_len615-seq_len)
}
dna615 <- cbind(dna615, info615)


### 4.2.2. Filter dataset ---------------------------------------------------

# Remove sequences with ambiguities or with less than 234 nt
rem <- dna615[dna615$N.ATCG < 615, "ID"]
dna615_red <- dna615[!(dna615$ID %in% rem),]

# Remove duplicated sequences
dup <- dna615_red[base::grepl("-", dna615_red$ID), "ID"]
dna615_red <- dna615_red[!(dna615_red$ID %in% dup), ]


### 4.1.3. Save clean fasta -------------------------------------------------

seq615_red <- dna615_red$Sequence
s615.red <- DNAStringSet(seq615_red)
names(s615.red) <- dna615_red$ID
writeXStringSet(s615.red, "fasta/Beluga_615bp_onlyATGC_n3106.fasta")


