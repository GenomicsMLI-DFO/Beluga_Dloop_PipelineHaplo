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
getwd()

# Clear workspace
rm(list = ls())

# Libraries
library(readxl)
library(dplyr)

# Functions
"%notin%" <- Negate("%in%")




# 1. Data -----------------------------------------------------------------

data234 <- read.csv("Sequences_Dloop234_n3314.csv")
str(data234)
data615 <- read.csv("Sequences_Dloop615_n3314.csv")
str(data615)
# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in 2a_HaploLibrary_234.R and 2b_HaploLibrary_615.R


# 2. Assign haplotype to each individual ----------------------------------

## 2.1. Upload haplotype libraries ------------------------------------------

lib234 <- read.csv('libraries/librairie_51_haplotypes234.csv')  # most recent haplotype library
colnames(lib234) <- c("hapl","seq")  # if it's not already the case
lib615 <- read.csv('libraries/librairie_137_haplotypes615.csv')  # most recent haplotype library
colnames(lib615) <- c("hapl","seq")  # if it's not already the case

## 2.2. Upload info on minimal sequence -----------------------------------

min_seq <- read.csv("polymorphismes_et_seq_minimale.csv", stringsAsFactors = F)  # table made in script 1
seq_start234 <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[2], " - "))[1])  # start of short minimal sequence
seq_stop234 <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[2], " - "))[2])  # end of short minimal sequence
seq_start615 <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[1], " - "))[1])  # start of long minimal sequence
seq_stop615 <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale[1], " - "))[2])  # end of long minimal sequence

## 2.3. Assign haplotype to each individual -------------------------------

hapind <- data.frame(matrix(ncol=2, nrow=0))
colnames(hapind) <- c("haplotype_234","haplotype_615")

seq234 <- toupper(data234$seq)
seq615 <- toupper(data615$seq)

for (i in 1:length(seq234)){
  if(data234$seq_utilisable[i] == 1) {
    if(substr(data234$seq[i], seq_start234, seq_stop234) %in% substr(lib234$seq, seq_start234, seq_stop234)) {
      hapind[i,1] <- lib234$hapl[which(substr(lib234$seq, seq_start234, seq_stop234) == substr(data234$seq[i], seq_start234, seq_stop234))]
    } else {
      # shows if haplotype is new relative to library: you MUST generate new haplo library if this is the case
      hapind[i,1] <- "unknown haplotype"
    }
  } else {
    # generates NAs for unusable sequences, to avoid confusion
    hapind[i,1] <- "NA"
  }
  if(data615$seq_utilisable[i] == 1) {
    if(substr(data615$seq[i], seq_start615, seq_stop615) %in% substr(lib615$seq, seq_start615, seq_stop615)) {
      hapind[i,2] <- lib615$hapl[which(substr(lib615$seq, seq_start615, seq_stop615) == substr(data615$seq[i], seq_start615, seq_stop615))]
    } else {
      # shows if haplotype is new relative to library: you MUST generate new haplo library if this is the case
      hapind[i,2] <- "unknown haplotype"
    }
  } else {
    # generates NAs for unusable sequences, to avoid confusion
    hapind[i,2] <- "NA"
  }
}

table(hapind$haplotype_234)
table(hapind$haplotype_615)
data <- merge(data234[,names(data234) %notin% c('seq','N.ATCG')], data615[,names(data615) %notin% c('seq','Numero_unique_extrait','N.ATCG')], by = "ID")
# table(data$Numero_unique_extrait.x == data$Numero_unique_extrait.y)
colnames(data) <- c('Numero_unique_specimen','Numero_unique_extrait','N_nucl_234','N_ambig_234','N_manquants_234','Sequence_utilisable_234',
                    'N_nucl_615','N_ambig_615','N_manquants_615','Sequence_utilisable_615')
data2 <- cbind(data, hapind)





# 3. Include new haplo in D-Loop sheet ------------------------------------

## 3.1. Upload ACCESS (D-Loop) dataset ------------------------------------

d <- read.csv("Dloop_MOBELS.csv")


## 3.2. Include new haplotypes in 'd' (D-Loop ACCESS) --------------------

data2$Numero_unique_specimen <- gsub("-.","",data2$Numero_unique_specimen)  # remove special identifier for duplicates

dloop <- left_join(d, data2, by = "Numero_unique_extrait")
table(dloop$Numero_unique_specimen.x == dloop$Numero_unique_specimen.y, useNA = "ifany")  # quality check: verify that specimen names are the same and there is no screw ups
# Keep Numero_unique_specimen.x when subsetting data frame since Numero_unique_specimen.y (data2) doesn't include all specimens

# Include info on library used to assign haplo to specimens - here since if done before the left join would introduce NA for specimens without sequence
dloop$Librairie_ref_234 <- paste('librairie',length(lib234$hapl),'haplotypes234',sep = '_')
dloop$Librairie_ref_615 <- paste('librairie',length(lib615$hapl),'haplotypes615',sep = '_')

# Subset dataset - same columns (and order) as original D-Loop sheet
dloop <- subset(dloop, select = c('Numero_unique_DLoop','Numero_unique_extrait','Numero_unique_specimen.x','Nom_Projet','Responsable_Dloop','Qualite_sequence',
                                  'Sequence_utilisable_234','Sequence_utilisable_615','No_run_F','No_plaque_F','No_puits_F','No_run_R','No_plaque_R','No_puits_R',
                                  'Sequence_consensus','N_nucl_234','N_ambig_234','N_manquants_234','haplotype_234','Librairie_ref_234','N_nucl_615','N_ambig_615',
                                  'N_manquants_615','haplotype_615','Librairie_ref_615','Modifications','Notes','Numero_unique_tissus','Numero_unique_extrait...22'))
colnames(dloop)[c(3)] <- c('Numero_unique_specimen')





# 4. Sequence quality: qualitative check ----------------------------------
# 0 unusable
# 1 usable without heterozygous nucleotides
# 2 usable with heterozygous nucleotides
# 3 incomplete
# 4 doubts on vaility of sequence
# 11 good sequence in other extraction

table(dloop$Qualite_sequence[dloop$Sequence_utilisable_234 %in% 0], useNA = 'ifany')
# 2    3   11
# 5    1    3
table(dloop$Qualite_sequence[data$Sequence_utilisable_234 %in% 1], useNA = 'ifany')
#   0    1    2    3    4   11 
# 206 3289   17    5   35   79

table(dloop$Qualite_sequence[dloop$Sequence_utilisable_615 %in% 0], useNA = 'ifany')
#  1  2  3 11 
# 19 13  4  5
table(dloop$Qualite_sequence[data$Sequence_utilisable_615 %in% 1], useNA = 'ifany')
#   0    1    2    3    4   11 
# 205 3257   17    5   33   78




# 5. Save dataset ---------------------------------------------------------

# Formatting some column values (only a few to repeat each time)
# str(dloop)
# table(dloop$No_run_F)
# dloop$No_run_F <- gsub("0.0", "0", dloop$No_run_F)
# table(dloop$No_plaque_F)
# dloop$No_plaque_F <- gsub("0.0", "0", dloop$No_plaque_F)
# table(dloop$No_puits_F)
# dloop$No_puits_F <- gsub("0.0", "0", dloop$No_puits_F)
# table(dloop$No_run_R)
# dloop$No_run_R <- gsub("0.0", "0", dloop$No_run_R)
# table(dloop$No_plaque_R)
# dloop$No_plaque_R <- gsub("0.0", "0", dloop$No_plaque_R)
# table(dloop$No_puits_R)
# dloop$No_puits_R <- gsub("0.0", "0", dloop$No_puits_R)
# table(dloop$Responsable_Dloop, useNA = "ifany")
# dloop[dloop$Responsable_Dloop %in% "0.0", "Responsable_Dloop"] <- NA
# dloop[dloop$Notes %in% "0.0", "Notes"] <- NA
# dloop[dloop$Responsable_Dloop %in% "0.0", "Responsable_Dloop"] <- NA
dloop[is.na(dloop$Sequence_utilisable_234), "Sequence_utilisable_234"] <- 0
dloop[is.na(dloop$Sequence_utilisable_615), "Sequence_utilisable_615"] <- 0

dloop <- arrange(dloop, Numero_unique_extrait)
write.csv(dloop, "Dloop_haplo_n3643.csv", row.names = F)  # Upload this directly on ACCESS file, D-Loop sheet
write.table(dloop, "Dloop_haplo_n3643.txt", sep = "\t" , row.names = F)
