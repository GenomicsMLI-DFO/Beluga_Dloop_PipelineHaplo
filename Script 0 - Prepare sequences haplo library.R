# Prepare sequences for updating haplotype libraries: long (HL, 615bp) and short (HS, 234 bp)
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

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")
library(Biostrings)
library(msa)

#library(ape)

# Data --------------------------------------------------------------------
# Multiple sequence alignment ---------------------------------------------

## Upload databases: originally in ACCESS folder on Drive
d <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20211216_MOBELS.xlsx", sheet = "D-Loop", na = "NA")  # database in ACCESS folder on Drive
s <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20211216_MOBELS.xlsx", sheet = "Specimens", na = "NA")

## Subset columns
# Dloop
#d <- data.frame(d)
str(d)
colnames(d)[3] <- "Numero_unique_extrait"
d <- subset(d, Qualite_sequence==1, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
d <- d[!is.na(d$Sequence_consensus),]
table(d$Qualite_sequence)

# Specimens
#s <- data.frame(s)
str(s)
s <- subset(s, select = c(Numero_unique_specimen, Nom_commun, Age, Sexe_visuel, Region_echantillonnage, Annee_echantillonnage:Jour_echantillonnage))

# Merge and remove narwhals and hybrids
dt <- merge(d, s, by = "Numero_unique_specimen")  # excludes specimens that are in s but not in d (18 specimens, nrow = 3641)
table(dt$Nom_commun)
dt <- subset(dt, Nom_commun == "Beluga", select = c(Numero_unique_specimen, Numero_unique_extrait, Sequence_consensus, N_nucl, Age:Jour_echantillonnage))

## Format Sequence_consensus
dt$Sequence_consensus <- toupper(dt$Sequence_consensus)
dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)
dt$N_nucl2 <- nchar(dt$Sequence_consensus)
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))
#     A      C      D      G      H      K      M      N      R      S      T      W      Y 
#630308 515581      1 336738      1     20     19      1     11     10 708704     41     21
# Don't remember what I wanted to do here... Interesting way the for loop to count for number of ambiguities in 1_Database...
#n <- data.frame(matrix(NA, nrow = nrow(dt), ncol = 16))
#colnames(n) <- c("Numero_unique_specimen","A","C","T","G","N","R","Y","K","M","S","W","B","D","H","V")
#for(i in 1:nrow(dt)){
#}

## Identify duplicated sequences
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
dloop$Numero_unique_specimen2 <- ifelse(is.na(dloop$Duplicated), dloop$Numero_unique_specimen,
                                        paste(dloop$Numero_unique_specimen, dloop$Duplicated, sep = "-"))
dloop <- subset(dloop, select = c(Numero_unique_specimen2, Numero_unique_extrait, Sequence_consensus, N_nucl, N_nucl2, Age:Jour_echantillonnage))
colnames(dloop)[1] <- "Numero_unique_specimen"
dloop <- arrange(dloop, Numero_unique_specimen)


# Multiple sequence alignment ---------------------------------------------
# Prepare DNAStringSet object
seq <- dloop$Sequence_consensus
#table(strsplit(seq[1], split = ""))  # verify which nucleotides or ambiguities are present (maybe do it for all sequence with a final summary table?)
dna <- DNAStringSet(seq, use.names = T)
names(dna) <- dloop$Numero_unique_specimen
#writeXStringSet(dna, "Beluga_n3284.fasta")

# Seq in revcomplement found in Mega
# S_20_00647 (duplicated), 01198, 01618, 01638, 02908, 03180, 03202
dna$S_20_00647 <- reverseComplement(dna$S_20_00647)
dna$S_20_01198 <- reverseComplement(dna$S_20_01198)
dna$S_20_01618 <- reverseComplement(dna$S_20_01618)
dna$S_20_01638 <- reverseComplement(dna$S_20_01638)
dna$S_20_02908 <- reverseComplement(dna$S_20_02908)
dna$S_20_03180 <- reverseComplement(dna$S_20_03180)
dna$S_20_03202 <- reverseComplement(dna$S_20_03202)
writeXStringSet(dna, "Beluga_complete_n3284.fasta")

# Alignment
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
writeXStringSet(alignment, "Beluga_alignment_complete_n3284.fasta")
#dna.algn <- readDNAStringSet("Beluga_alignment_complete_n3284.fasta")  # upload complete alignment

# Subset sequences
# Use bp strings at start and end of consensus sequence to define the cutting loci.
# Define F and R 'primers
cs <- readDNAStringSet(filepath = "~/Documents/Post-Docs/IML/MOBELS/dloop/Fasta/Ref sequence/Sequence_Dloop_Ref_615pb_consensus.fasta")
cs
dloop.F <- "ACTACGTCAGTATTAAATAAA"
dloop.R <- "GCTGGACCTGTGTGTATTTTT"

# Define cutting position
res.dloop.F <- vmatchPattern(DNAString(dloop.F), DNAStringSet(dna.algn), max.mismatch = 4)
cut.dloop.F.int <- res.dloop.F@ends %>% unlist() %>% table() %>% dimnames(.)
cut.dloop.F <- as.numeric(as.character(cut.dloop.F.int[[1]])) - nchar(dloop.F) + 1

res.dloop.R <- vmatchPattern(DNAString(dloop.R), DNAStringSet(dna.algn), max.mismatch = 4)
cut.dloop.R.int <- res.dloop.R@ends %>% unlist() %>% table() %>% dimnames(.)
cut.dloop.R <- as.numeric(as.character(cut.dloop.R.int[[1]]))

Dloop <- subseq(DNAStringSet(dna.algn), start = cut.dloop.F, end = cut.dloop.R)
print(Dloop, show = "complete")
table(Dloop@ranges@width)
Dloop
s615 <- Dloop
names(s615) <- names(Dloop)
writeXStringSet(s615, "Beluga_615bp_n3284.fasta")
#s615 <- readDNAStringSet("Beluga_615bp_n3284.fasta")  # upload 615bp alignment


# Add metadata to DNAStringSet and save dataset ---------------------------
#metadata(DNA.dloop)$Meta_region <- NA
dloop.dna <- data.frame(ID = names(Dloop),
                        Sequence = Dloop)
dloop.dna <- left_join(dloop.dna, dloop[,c("Numero_unique_specimen","Numero_unique_extrait","Age","Sexe_visuel","Region_echantillonnage",
                                           "Annee_echantillonnage","Mois_echantillonnage","Jour_echantillonnage")], by = c("ID"="Numero_unique_specimen"))
write.table(dloop.dna, file = "Sequences_Dloop_all_n3284.txt", row.names = F)


# Prepare dataset for script 1 --------------------------------------------
# Clean sequences needed: 615 bp and NO ambiguous nt
dloop.dna
seq <- dloop.dna$Sequence
nt <- c("A","T","C","G")
ambiguous <- c("N","R","Y","K","M","S","W","B","D","H","V")
exp_seq_len <- 615

info <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq)){
    if(is.na(seq[i])){
        info[i,] <- c("NA","NA","NA","NA")	
    }else{
        seq_len <- sum(str_count(seq[i],c(nt, ambiguous)))
    }	
    info[i,] <- c(seq_len, sum(str_count(seq[i], nt)), sum(str_count(seq[i], ambiguous)), exp_seq_len-seq_len)
}
dloop.dna <- cbind(dloop.dna, info)

rem <- dloop.dna[dloop.dna$N.ATCG < 615, "ID"]
dloop.dna_red <- dloop.dna[!(dloop.dna$ID %in% rem),]

# Remove duplicated sequences
dup <- dloop.dna_red[base::grepl("-", dloop.dna_red$ID), "ID"]
dloop.dna_red <- dloop.dna_red[!(dloop.dna_red$ID %in% dup), ]

# Save clean fasta
seq_red <- dloop.dna_red$Sequence
s615.red <- DNAStringSet(seq_red)
names(s615.red) <- dloop.dna_red$ID
writeXStringSet(s615.red, "Beluga_615bp_onlyATGC_n3102.fasta")





### Complete dloop seuquence: https://www.ncbi.nlm.nih.gov/nuccore/U18117.1
# Cut something before 129 and somewhat later 355 (SNPs position found in Brown Gladden et al. 1997)
## 226 bp between 129 and 355, the sequence in Brown Gladden et al. 1997 is 234 bp
