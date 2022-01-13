# Info --------------------------------------------------------------------

# Prepare sequences to update haplotype libraries: long (HL, 615bp) and short (HS, 234 bp)
# Multiple sequences alignment
# 
# 
# Luca Montana
# 2021-12-17
#

getwd()
rm(list = ls())


# Libraries and function --------------------------------------------------
if(!require(tidyverse)){install.packages("tidyverse")}
library(readxl)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("msa")
library(Biostrings)
library(msa)



# Data --------------------------------------------------------------------

## Upload databases: originally in ACCESS folder on Drive -----------------
d <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20220112_MOBELS.xlsx", sheet = "D-Loop", na = "NA")  # database in ACCESS folder on Drive
s <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20220112_MOBELS.xlsx", sheet = "Specimens", na = "NA")

## Format input database for MSA ------------------------------------------
# Dloop database
str(d)  # 3643 rows
colnames(d)[2] <- "Numero_unique_extrait"
d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
#d <- subset(d, Qualite_sequence==1, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))  # removes 341 rows
d <- d[!is.na(d$Sequence_consensus),]  # removes 1 additional sequence (315 rows if not removing sequences quality !=1)
table(d$Qualite_sequence)
#   1    2    3    4   11 
#3300   17    5    1    5
table(is.na(d$Qualite_sequence))  # No NAs in Qualite_sequence

# Specimens databse
str(s)
colnames(s)[3] <- "Numero_unique_specimen"
s <- s[,c("Numero_unique_specimen","Nom_commun")]

# Merge and remove narwhals and hybrids
dt <- merge(d, s, by = "Numero_unique_specimen")
table(dt$Nom_commun)
dt <- dt[dt$Nom_commun %in% "Beluga", c("Numero_unique_specimen","Numero_unique_extrait","Sequence_consensus","N_nucl")]  # removes 14 specimens

### Format Sequence_consensus ---------------------------------------------
dt$Sequence_consensus <- toupper(dt$Sequence_consensus)
dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)
dt$N_nucl2 <- nchar(dt$Sequence_consensus)
table(dt$N_nucl2)
dt <- dt[!dt$N_nucl2 < 200,]  # remove any sequence that is very short (there should be none)

# Verify which ambiguities are included in sequences at present
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))


### Identify duplicated sequences -----------------------------------------
# Necessary to recognise who's who in fasta files
dt$Duplicated <- NA
uni <- dt[!duplicated(dt$Numero_unique_specimen),]
dup <- dt[duplicated(dt$Numero_unique_specimen),]
dup2 <- dup[duplicated(dup$Numero_unique_specimen),]
dup <- dup[!duplicated(dup$Numero_unique_specimen),]
dup$Duplicated <- 2
dup3 <- dup2[duplicated(dup2$Numero_unique_specimen),]
dup2 <- dup2[!duplicated(dup2$Numero_unique_specimen),]
dup2$Duplicated <- 3
dup3$Duplicated <- 4
dloop <- rbind(uni,dup,dup2,dup3)

# Add identifier in Numero_unique_specimen
dloop$Numero_unique_specimen2 <- ifelse(is.na(dloop$Duplicated), dloop$Numero_unique_specimen,
                                        paste(dloop$Numero_unique_specimen, dloop$Duplicated, sep = "-"))
dloop <- subset(dloop, select = c(Numero_unique_specimen2, Numero_unique_extrait, Sequence_consensus, N_nucl, N_nucl2))
colnames(dloop)[1] <- "Numero_unique_specimen"
dloop <- arrange(dloop, Numero_unique_specimen)


# Multiple sequence alignment ---------------------------------------------

## Prepare DNAStringSet object --------------------------------------------
seq <- dloop$Sequence_consensus
#table(strsplit(seq[1], split = ""))  # verify which nucleotides or ambiguities are present (maybe do it for all sequence with a final summary table?)
dna <- DNAStringSet(seq)
names(dna) <- dloop$Numero_unique_specimen
#writeXStringSet(dna, "Beluga_n3283.fasta")

## Seq in revcomplement found in Mega -------------------------------------
# S_20_00647 (duplicated), 01198, 01618, 01638, 02908, 03180, 03202
dna$S_20_00647 <- reverseComplement(dna$S_20_00647)
dna$S_20_01198 <- reverseComplement(dna$S_20_01198)
dna$S_20_01618 <- reverseComplement(dna$S_20_01618)
dna$S_20_01638 <- reverseComplement(dna$S_20_01638)
dna$S_20_02908 <- reverseComplement(dna$S_20_02908)
dna$S_20_03180 <- reverseComplement(dna$S_20_03180)
dna$S_20_03202 <- reverseComplement(dna$S_20_03202)
writeXStringSet(dna, "Beluga_complete_seq_n3314.fasta")

## MSA --------------------------------------------------------------------
dna.algn <- msa(dna, method = "Muscle", gapOpening = 10000, gapExten = 400, maxiters = 30, type = "dna",
                order = "input", verbose = T)
print(dna.algn, show = "complete")  # 87 to 701+1
#checkAlignment(dna.algn)  # ape function, need a DNAbin object
# msaPrettyPrint(dna.algn, output = "pdf", file = "beluga_all_msa_1-58.pdf", subset=c(1:58), showConsensus = "none", showNames = "left")
# msaPrettyPrint(dna.algn, output = "pdf", file = "beluga_all_msa_59-116.pdf", subset=c(59:116), showConsensus = "none", showNames = "left")
# msaPrettyPrint(dna.algn, output = "pdf", file = "beluga_all_msa_117-174.pdf", subset=c(117:174), showConsensus = "none", showNames = "left")
# MegaX (remember to add consensus sequence 615 bp at the top)
# Settings: Gap open = -10000.00; Gap Extend = -400.00; Max Iterations = 30; Cluster Method = UPGMA
alignment <- DNAStringSet(dna.algn)  # to save the alignment
writeXStringSet(alignment, "Beluga_alignment_complete_n3314.fasta")
#dna.algn <- readDNAStringSet("Beluga_alignment_complete_n3314.fasta")  # upload complete alignment


# Cut sequences -----------------------------------------------------------
# Use bp strings at start and end of consensus sequences (615bp and 234bp) to define the cutting loci

## Cut sequences 615 bp ---------------------------------------------------

### Define F and R 'primers' ----------------------------------------------
cr615 <- readDNAStringSet(filepath = "~/Documents/Post-Docs/IML/MOBELS/dloop/Fasta/Ref sequence/Sequence_Dloop_Ref_615pb_consensus.fasta")
F615 <- toString(subseq(cr615, start = 1, width = 21))  # at the beginning of 5'-end of CR (nt = 38)
R615 <- toString(subseq(cr615, end = 615, width = 21))


### Define cutting position -----------------------------------------------
res.F615 <- vmatchPattern(DNAString(F615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F615.int <- res.F615@ends %>%
    unlist() %>%
    table() %>%
    dimnames(.)
cut.F615 <- as.integer(as.character(cut.F615.int[[1]])) - nchar(F615) + 1

res.R615 <- vmatchPattern(DNAString(R615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R615.int <- res.R615@ends %>%
    unlist() %>%
    table() %>%
    dimnames(.)
cut.R615 <- as.numeric(as.character(cut.R615.int[[1]]))

### Cut sequences to 615 bp and save fasta ---------------------------------
Dloop615 <- subseq(DNAStringSet(dna.algn), start = cut.F615, end = cut.R615)
print(Dloop615, show = "complete")
table(Dloop615@ranges@width)
Dloop615
writeXStringSet(Dloop615, "Beluga_615bp_n3314.fasta")
#Dloop615 <- readDNAStringSet("Beluga_615bp_n3314.fasta")  # upload 615bp alignment

### Save dataset ------------------------------------------------------------
#metadata(DNA.dloop)$Meta_region <- NA
dna615 <- data.frame(ID = names(Dloop615),
                     Sequence = Dloop615)
dna615 <- left_join(dna615, dloop[,c("Numero_unique_specimen","Numero_unique_extrait")], by = c("ID"="Numero_unique_specimen"))
write.table(dna615, file = "Sequences_Dloop615_n3314.txt", row.names = F)


## Cut sequences 234 bp ----------------------------------------------------

### Define F and R 'primers' -----------------------------------------------
# In the 235 basepairs routinely sequenced (position 134 to 384) we have found 18 variable sites
# expressing 37 haplotypes in 450 individuals. FROM Lillie et al. 1996.
# Nineteen polymorphic nucleotide positions were found and nearly all changes were transitions
# [first SNP = 129; last SNP = 355]. FROM Brown Gladden et al. 1997.
# The mtDNA locus we used consists of 234 nucleotides that are found at the beginning of the d-loop
# region of the molecule (Brown Gladden et al. 1997). FROM De March & Postma 2003.
# Primer 5'-3' used in Brown Gladden et al. and De March & Postma:
# Bel5 <- "ACATTTTACTGTGACTATTG"  # at the beginning of 5'-end of CR (nt = 71)
cr917 <- readDNAStringSet(filepath = "~/Documents/Post-Docs/IML/MOBELS/dloop/Fasta/Ref sequence/Sequence_Dloop_complete.fasta")
F234 <- toString(subseq(cr917, start = 126, width = 21))  # at the beginning of 5'-end of CR (nt = 38)
R234 <- toString(subseq(cr917, end = 359, width = 21))


### Define cutting position ------------------------------------------------
res.F234 <- vmatchPattern(DNAString(F234), DNAStringSet(dna.algn), max.mismatch = 2)
cut.F234.int <- res.F234@ends %>%
    unlist() %>%
    table() %>%
    dimnames(.)
cut.F234 <- as.numeric(as.character(cut.F234.int[[1]])) - nchar(F234) + 1 # nt 915: normal as alignment starts BEFORE the 5'-strand of CR

res.R234 <- vmatchPattern(DNAString(R234), DNAStringSet(dna.algn), max.mismatch = 4)
cut.R234.int <- res.R234@ends %>%
    unlist() %>%
    table() %>%
    dimnames(.)
cut.R234 <- as.numeric(as.character(cut.R234.int[[1]]))

### Cut sequences to 234 bp and save fasta ----------------------------------
Dloop234 <- subseq(DNAStringSet(dna.algn), start = cut.F234, end = cut.R234)
print(Dloop234, show = "complete")
table(Dloop234@ranges@width)
writeXStringSet(Dloop234, "Beluga_234bp_n3314.fasta")
#Dloop234 <- readDNAStringSet("Beluga_234bp_n3314.fasta")  # upload 234bp alignment

### Save dataset ------------------------------------------------------------
#metadata(DNA.dloop)$Meta_region <- NA
dna234 <- data.frame(ID = names(Dloop234),
                     Sequence = Dloop234)
dna234 <- left_join(dna234, dloop[,c("Numero_unique_specimen","Numero_unique_extrait")], by = c("ID"="Numero_unique_specimen"))
write.table(dna234, file = "Sequences_Dloop234_n3314.txt", row.names = F)


# Prepare datasets (615bp and 234bp) for script 1 ---------------------------
# Clean sequences needed: 615/234 bp and NO ambiguous nt
nt <- c("A","T","C","G")
ambiguous <- c("N","R","Y","K","M","S","W","B","D","H","V")

## 615 bp -------------------------------------------------------------------

### Include no nt, no ATCG, no ambig and no missing nt in dataset -----------
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

### Remove specimens with number nt < 615 -----------------------------------
rem <- dna615[dna615$N.ATCG < 615, "ID"]
dna615_red <- dna615[!(dna615$ID %in% rem),]

### Remove duplicated sequences ---------------------------------------------
dup <- dna615_red[base::grepl("-", dna615_red$ID), "ID"]
dna615_red <- dna615_red[!(dna615_red$ID %in% dup), ]

### Save clean fasta 615 ----------------------------------------------------
seq615_red <- dna615_red$Sequence
s615.red <- DNAStringSet(seq615_red)
names(s615.red) <- dna615_red$ID
writeXStringSet(s615.red, "Beluga_615bp_onlyATGC_n3107.fasta")

## 234 bp -------------------------------------------------------------------

### Include no nt, no ATCG, no ambig and no missing nt in dataset -----------
seq234 <- dna234$Sequence
exp_seq_len234 <- 234

info234 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info234) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq234)){
    if(is.na(seq234[i])){
        info234[i,] <- c("NA","NA","NA","NA")	
    }else{
        seq_len <- sum(str_count(seq234[i],c(nt, ambiguous)))
    }	
    info234[i,] <- c(seq_len, sum(str_count(seq234[i], nt)), sum(str_count(seq234[i], ambiguous)), exp_seq_len234-seq_len)
}
dna234 <- cbind(dna234, info234)

### Remove specimens with number nt < 234 -----------------------------------
rem <- dna234[dna234$N.ATCG < 234, "ID"]  # none
dna234_red <- dna234[!(dna234$ID %in% rem),]

### Remove duplicated sequences ---------------------------------------------
dup <- dna234_red[base::grepl("-", dna234_red$ID), "ID"]
dna234_red <- dna234_red[!(dna234_red$ID %in% dup), ]

### Save clean fasta 234 ----------------------------------------------------
seq234_red <- dna234_red$Sequence
s234.red <- DNAStringSet(seq234_red)
names(s234.red) <- dna234_red$ID
writeXStringSet(s234.red, "Beluga_234bp_onlyATGC_n3176.fasta")

