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



# Data --------------------------------------------------------------------

## Upload databases: originally in ACCESS folder on Drive -----------------
d <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20211216_MOBELS.xlsx", sheet = "D-Loop", na = "NA")  # database in ACCESS folder on Drive
s <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20211216_MOBELS.xlsx", sheet = "Specimens", na = "NA")

## Format input database for MSA ------------------------------------------
# Dloop database
str(d)
colnames(d)[3] <- "Numero_unique_extrait"
d <- subset(d, Qualite_sequence==1, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
d <- d[!is.na(d$Sequence_consensus),]
table(d$Qualite_sequence)

# Specimens databse
str(s)
s <- s[,c("Numero_unique_specimen","Nom_commun","Age","Sexe_visuel","Region_echantillonnage","Annee_echantillonnage","Mois_echantillonnage","Jour_echantillonnage")]

# Merge and remove narwhals and hybrids
dt <- merge(d, s, by = "Numero_unique_specimen")  # excludes specimens that are in s but not in d (18 specimens, nrow = 3641)
table(dt$Nom_commun)
dt <- dt[dt$Nom_commun %in% "Beluga", c("Numero_unique_specimen","Numero_unique_extrait","Sequence_consensus","N_nucl","Age","Sexe_visuel",
                                        "Region_echantillonnage","Annee_echantillonnage","Mois_echantillonnage","Jour_echantillonnage")]

### Format Sequence_consensus ---------------------------------------------
dt$Sequence_consensus <- toupper(dt$Sequence_consensus)
dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)
dt$N_nucl2 <- nchar(dt$Sequence_consensus)

# Verify which ambiguities are included in sequences at present
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))
#     A      C      D      G      H      K      M      N      R      S      T      W      Y 
#630308 515581      1 336738      1     20     19      1     11     10 708704     41     21

### Identify duplicated sequences -----------------------------------------
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
dloop <- subset(dloop, select = c(Numero_unique_specimen2, Numero_unique_extrait, Sequence_consensus, N_nucl, N_nucl2, Age:Jour_echantillonnage))
colnames(dloop)[1] <- "Numero_unique_specimen"
dloop <- arrange(dloop, Numero_unique_specimen)


# Multiple sequence alignment ---------------------------------------------

## Prepare DNAStringSet object --------------------------------------------
seq <- dloop$Sequence_consensus
#table(strsplit(seq[1], split = ""))  # verify which nucleotides or ambiguities are present (maybe do it for all sequence with a final summary table?)
dna <- DNAStringSet(seq)
names(dna) <- dloop$Numero_unique_specimen
#writeXStringSet(dna, "Beluga_n3284.fasta")

## Seq in revcomplement found in Mega -------------------------------------
# S_20_00647 (duplicated), 01198, 01618, 01638, 02908, 03180, 03202
dna$S_20_00647 <- reverseComplement(dna$S_20_00647)
dna$S_20_01198 <- reverseComplement(dna$S_20_01198)
dna$S_20_01618 <- reverseComplement(dna$S_20_01618)
dna$S_20_01638 <- reverseComplement(dna$S_20_01638)
dna$S_20_02908 <- reverseComplement(dna$S_20_02908)
dna$S_20_03180 <- reverseComplement(dna$S_20_03180)
dna$S_20_03202 <- reverseComplement(dna$S_20_03202)
writeXStringSet(dna, "Beluga_complete_n3284.fasta")

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
writeXStringSet(alignment, "Beluga_alignment_complete_n3284.fasta")
#dna.algn <- readDNAStringSet("Beluga_alignment_complete_n3284.fasta")  # upload complete alignment


# Cut sequences -----------------------------------------------------------
# Use bp strings at start and end of consensus sequences (615bp and 234bp) to define the cutting loci

## Cut sequences 615 bp ---------------------------------------------------

### Define F and R 'primers' ----------------------------------------------
cr615 <- readDNAStringSet(filepath = "~/Documents/Post-Docs/IML/MOBELS/dloop/Fasta/Ref sequence/Sequence_Dloop_Ref_615pb_consensus.fasta")
cr615
F615 <- "ACTACGTCAGTATTAAATAAA"  # at the beginning of 5'-end of CR (nt = 38)
R615 <- "GCTGGACCTGTGTGTATTTTT"

### Define cutting position -----------------------------------------------
res.F615 <- vmatchPattern(DNAString(F615), DNAStringSet(dna.algn), max.mismatch = 4)
cut.F615.int <- res.F615@ends %>%
    unlist() %>%
    table() %>%
    dimnames(.)
cut.F615 <- as.numeric(as.character(cut.F615.int[[1]])) - nchar(F615) + 1

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
writeXStringSet(Dloop615, "Beluga_615bp_n3284.fasta")
#Dloop615 <- readDNAStringSet("Beluga_615bp_n3284.fasta")  # upload 615bp alignment

### Save dataset ------------------------------------------------------------
#metadata(DNA.dloop)$Meta_region <- NA
dna615 <- data.frame(ID = names(Dloop615),
                     Seq615 = Dloop615)
dna615 <- left_join(dna615, dloop[,c("Numero_unique_specimen","Numero_unique_extrait","Age","Sexe_visuel","Region_echantillonnage",
                                     "Annee_echantillonnage","Mois_echantillonnage","Jour_echantillonnage")], by = c("ID"="Numero_unique_specimen"))
write.table(dna615, file = "Sequences_Dloop615_all_n3284.txt", row.names = F)


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
cr917 <- as.character(cr917$`U18117.1 Delphinapterus leucas mitochondrion control region, complete sequence`)
F234 <- substring(cr917, first = 126, last = 146)  # A few nt before the first SNP found by Brown Gladden et al . 1997
R234 <- substring(cr917, first = 339, last = 359)  # A few nt after the last SNP found by Brown Gladden et al . 1997

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
writeXStringSet(Dloop234, "Beluga_234bp_n3284.fasta")
#Dloop234 <- readDNAStringSet("Beluga_234bp_n3284.fasta")  # upload 615bp alignment

### Save dataset ------------------------------------------------------------
#metadata(DNA.dloop)$Meta_region <- NA
dna234 <- data.frame(ID = names(Dloop234),
                     Seq234 = Dloop234)
dna234 <- left_join(dna234, dloop[,c("Numero_unique_specimen","Numero_unique_extrait","Age","Sexe_visuel","Region_echantillonnage",
                                     "Annee_echantillonnage","Mois_echantillonnage","Jour_echantillonnage")], by = c("ID"="Numero_unique_specimen"))
write.table(dna234, file = "Sequences_Dloop234_all_n3284.txt", row.names = F)


# Prepare datasets (615bp and 234bp) for script 1 --------------------------------------------
# Clean sequences needed: 615 bp and NO ambiguous nt
seq615 <- as.character(Dloop615)
seq234 <- as.character(Dloop234)
nt <- c("A","T","C","G")
ambiguous <- c("N","R","Y","K","M","S","W","B","D","H","V")
exp_seq_len615 <- 615
exp_seq_len234 <- 234

info615 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(info) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(seq615)){
    if(is.na(seq615[i])){
        info[i,] <- c("NA","NA","NA","NA")	
    }else{
        seq_len <- sum(str_count(seq615[i],c(nt, ambiguous)))
    }	
    info[i,] <- c(seq_len, sum(str_count(seq615[i], nt)), sum(str_count(seq615[i], ambiguous)), exp_seq_len615-seq_len)
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









