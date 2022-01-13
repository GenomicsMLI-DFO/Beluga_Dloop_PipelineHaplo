# Prepare sequence dataset for haplotype assignment:
# Compiling the d_haplotypes library
# 
# Benjamin Hornoy
#
#

getwd()
rm(list = ls())

# Libraries and function --------------------------------------------------
library(readxl)


# Data --------------------------------------------------------------------
data <- read.csv("Sequences_Dloop234_n3314.csv")
str(data)
# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in Script 2 - Compiling d_haplotype library.R


# Assign haplotype to each individuals ------------------------------------
lib <- read.csv('librairie_51_haplotypes234.csv') # most recent haplotype library
colnames(lib) <- c("hapl","seq")  # if it's not already the case

## Upload info on minimal sequence ----------------------------------------
min_seq <- read.csv("polymorphismes_et_seq_minimale.csv", stringsAsFactors = F)  # table made in script 1
seq_start <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[2], " - "))[1])  # start of minimal sequence
seq_stop <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[2], " - "))[2])  # end of minimal sequence

## Upload info Qualite_sequence -------------------------------------------
# 0 unusable
# 1 usable without heterozygous nucleotides
# 2 usable with heterozygous nucleotides
# 3 incomplete
# 4 doubts on vaility of sequence
# 11 good sequence in other extraction
d <- read_excel("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/ACCESS/20220112_MOBELS.xlsx", sheet = "D-Loop", na = "NA")
colnames(d)[2] <- "Numero_unique_extrait"
data <- merge(data, d[,c("Numero_unique_extrait","Qualite_sequence")], by = "Numero_unique_extrait")
table(data$Qualite_sequence[data$seq_utilisable=="no"])
#2  3 11 
#5  1  3
table(data$Qualite_sequence[data$seq_utilisable=="yes"])
#   1    2    3    4   11 
#3286   12    4    1    2 

## Assign haplotype to each individual ------------------------------------
hapind <- data.frame(matrix(ncol=1, nrow=0))
colnames(hapind) <- c("haplotype")

sequences <- toupper(data$seq)
for (i in 1:length(sequences)){
    if(data$seq_utilisable[i] == "yes") {
        if(substr(data$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {
            hapind[i,1] <- lib$hapl[which(substr(lib$seq, seq_start, seq_stop) == substr(data$seq[i], seq_start, seq_stop))]
        } else {
            # shows if haplotype is new relative to library: you MUST generate new haplo library if this is the case
            hapind[i,1] <- "unknown haplotype"
        }
    } else {
        # generates NAs for unusable sequences, to avoid confusion
        hapind[i,1] <- "NA"
    }
}

table(hapind)
data2 <- cbind(data, hapind)

write.csv(data2, "Dloop_haplo234_n3314.csv", row.names=F)

