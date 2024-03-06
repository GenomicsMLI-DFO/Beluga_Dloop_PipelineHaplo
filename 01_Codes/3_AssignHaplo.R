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

data234 <- read.csv("00_Data/02_dloop_clean/Sequences_Dloop234_n3441.csv")
str(data234)
data615 <- read.csv("00_Data/02_dloop_clean/Sequences_Dloop615_n3441.csv")
str(data615)

# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in 2a_HaploLibrary_234.R and 2b_HaploLibrary_615.R

# 2. Assign haplotype to each individual ----------------------------------

## 2.1. Upload haplotype libraries ------------------------------------------

lib234 <- read.csv('02_Results/00_libraries/librairie_54_haplotypes234.csv')  # most recent haplotype library
colnames(lib234) <- c("hapl","seq")  # if it's not already the case
lib615 <- read.csv('02_Results/00_libraries/librairie_143_haplotypes615.csv')  # most recent haplotype library
colnames(lib615) <- c("hapl","seq")  # if it's not already the case

## 2.2. Upload info on minimal sequence -----------------------------------

min_seq <- read.csv("02_Results/01_poly_seq_min/polymorphismes_et_seq_minimale.csv", stringsAsFactors = F)  # table made in script 1
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
    hapind[i,1] <- NA
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
    hapind[i,2] <- NA
  }
}
# if two sequences have same minimal sequence it's a consequence of the update of the haplo library were a false SNP has been corrected and minimal sequence has shrank
# as a consequence

table(hapind$haplotype_234)
table(hapind$haplotype_615)
data <- merge(data234[,names(data234) %nin% c('seq','N.ATCG')], data615[,names(data615) %nin% c('seq','Numero_unique_extrait','N.ATCG')], by = "ID")
data <- subset(data, select = c(ID,Numero_unique_extrait,No_plaque_sequencage_F.x,No_puits_sequencage_F.x,No_plaque_sequencage_R.x,No_puits_sequencage_R.x,
                                N.nucl.x,N.ambig.x,N.manquants.x,seq_utilisable.x,N.nucl.y,N.ambig.y,N.manquants.y,seq_utilisable.y))
colnames(data) <- c('Numero_unique_specimen','Numero_unique_extrait','No_plaque_sequencage_F','No_puits_sequencage_F','No_plaque_sequencage_R',
                    'No_puits_sequencage_R','N_nucl_234','N_ambig_234','N_manquants_234','Sequence_utilisable_234','N_nucl_615','N_ambig_615',
                    'N_manquants_615','Sequence_utilisable_615')
data2 <- cbind(data, hapind)
data2$Numero_unique_specimen <- gsub("-.","",data2$Numero_unique_specimen)  # remove special identifier for duplicates

# Some specimens are duplicated and the same DNA extraction is the source: possible issue later when creating dloop dataframe.
# Adding more metadata (plate and well numbers) avoids duplication
table(duplicated(data2[,c('Numero_unique_specimen','Numero_unique_extrait','No_plaque_sequencage_F','No_puits_sequencage_F',
                          'No_plaque_sequencage_R','No_puits_sequencage_R')]))
# dup_spec <- data2$Numero_unique_specimen[duplicated(data2[,c('Numero_unique_specimen','Numero_unique_extrait')])]
# dup_ext <- data2$Numero_unique_extrait[duplicated(data2[,c('Numero_unique_specimen','Numero_unique_extrait')])]
# dup <- data2[data2$Numero_unique_specimen %in% dup_spec & data2$Numero_unique_extrait %in% dup_ext,]




# 3. Include new haplo in D-Loop sheet ------------------------------------

## 3.1. Upload ACCESS (D-Loop) dataset ------------------------------------

d <- read_excel("../ACCESS/20221202_Mobels_Chlamys.xlsx", sheet = "Sequencage", na = "NA",    # remember to specify right path to beluga ACCESS dataset
                col_types = c(rep("text",13),rep("numeric",3),rep("text",2),rep("numeric",1),rep("text",3)))
# colnames(d)[2] <- "Numero_unique_extrait"

## 3.2. Include new haplotypes in 'd' (D-Loop ACCESS) --------------------

dloop <- left_join(d, data2, by = c('Numero_unique_extrait','No_plaque_sequencage_F','No_puits_sequencage_F','No_plaque_sequencage_R','No_puits_sequencage_R'))
table(dloop$Numero_unique_specimen.x == dloop$Numero_unique_specimen.y, useNA = "ifany")  # quality check: verify that specimen names are the same and there is no screw ups
# Keep Numero_unique_specimen.x when subsetting data frame since Numero_unique_specimen.y (data2) doesn't include all specimens

# Include info on library used to assign haplo to specimens - here since if done before the left join would introduce NA for specimens without sequence
dloop$Librairie_ref_234 <- paste('librairie',length(lib234$hapl),'haplotypes234',sep = '_')
dloop$Librairie_ref_615 <- paste('librairie',length(lib615$hapl),'haplotypes615',sep = '_')

# Subset dataset - same columns (and order) as original D-Loop sheet
dloop <- subset(dloop, select = c('Numero_unique_Dloop','Numero_unique_extrait','Numero_unique_specimen.x','Nom_Projet','Responsable_Dloop','No_plaque_F',
                                  'No_puits_F','No_plaque_R','No_puits_R','Sequence_consensus','N_nucl_615','N_ambig_615','N_manquants_615','haplotype_615',
                                  'Librairie_ref_615','Sequence_utilisable_615','Modifications_Dloop','Notes_Dloop','Numero_unique_tissus','Numero_unique_extrait_2'))
colnames(dloop)[c(3,11:16)] <- c('Numero_unique_specimen','N_nucl','N_ambig','N_manquants','Haplotype','Librairie_ref','Used_for_librairie_ref')
# dloop <- subset(dloop, select = c('Numero_unique_Dloop','Numero_unique_extrait','Numero_unique_specimen.x','Nom_Projet','Responsable_Dloop','No_plaque_F',
#                                   'No_puits_F','No_plaque_R','No_puits_R','Sequence_consensus','N_nucl_234','N_ambig_234','N_manquants_234','haplotype_234',
#                                   'Librairie_ref_234','N_nucl_615','N_ambig_615','N_manquants_615','haplotype_615','Librairie_ref_615','Modifications','Notes',
#                                   'Numero_unique_tissus','Numero_unique_extrait_2'))
# colnames(dloop)[c(3)] <- c('Numero_unique_specimen')
dloop$Used_for_librairie_ref[is.na(dloop$Used_for_librairie_ref)] <- 0




# # 4. Sequence quality: qualitative check ----------------------------------
# # 0 unusable
# # 1 usable without heterozygous nucleotides
# # 2 usable with heterozygous nucleotides
# # 3 incomplete
# # 4 doubts on vaility of sequence
# # 11 good sequence in other extraction
# 
# table(dloop$Qualite_sequence[dloop$Sequence_utilisable_234 %in% 0], useNA = 'ifany')
# # 2    3   11
# # 5    1    3
# table(dloop$Qualite_sequence[data$Sequence_utilisable_234 %in% 1], useNA = 'ifany')
# #   0    1    2    3    4   11
# # 206 3289   17    5   35   79
# 
# table(dloop$Qualite_sequence[dloop$Sequence_utilisable_615 %in% 0], useNA = 'ifany')
# #  1  2  3 11
# # 19 13  4  5
# table(dloop$Qualite_sequence[data$Sequence_utilisable_615 %in% 1], useNA = 'ifany')
# #   0    1    2    3    4   11
# # 205 3257   17    5   33   78




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
# dloop[is.na(dloop$Sequence_utilisable_234), "Sequence_utilisable_234"] <- 0
# dloop[is.na(dloop$Sequence_utilisable_615), "Sequence_utilisable_615"] <- 0
dloop$Sequence_consensus <- gsub("-", "", dloop$Sequence_consensus)  # No breaks within sequences, remove '-' on the edges. Necessary step to avoid losing part of sequences
# that starts with "---". Apparently, excel and csv file are cut after about 255 characters if they start with "---"

dloop <- arrange(dloop, Numero_unique_extrait)
write.csv(dloop, "02_Results/02_ACCESS/Dloop_haplo_n3688.csv", row.names = F)  # Once exported, save this file as .xlsx, then copy-paste its content in a new sheet that will take the place
# of D-Loop shhet on MOBELS ACCESS file. Once this is done, run a few qualitative checks (notably on sequences lengths) to see that the export was done correctly (see
# LL 181-182 for an explanation - for example S_20_02438 sequence starts with "--GATT" to align it to other sequences. This sequence length is of 716 characters if 
# including "--" or 697 when excluding "--" as there are some at the end of the sequence as well. However, if saved with "--" excel cuts it to about 255 characters).



