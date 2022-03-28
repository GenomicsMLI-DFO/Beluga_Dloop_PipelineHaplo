# Info --------------------------------------------------------------------
# 
# Authors: Benjamin Hornoy, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-12-17
# 
# Overview: Bonus script to remove spaces and punctuation from D-Loop sheet (Beluga dataset).
# It should not be needed anymore - unless new metadata come in and have unfriendly format for R, ACCESS or conversion csv-to-excel
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
# library(dplyr)
# library(data.table)  # rleid function

# Functions
# "%notin%" <- Negate("%in%")




# 1. Upload data ----------------------------------------------------------

# Originally in ACCESS folder on Drive. Specify the path to the directory where the file is stored
d <- read_excel("../ACCESS/20220318_MOBELS.xlsx", sheet = "D-Loop", na = "NA")
str(d)
colnames(d)[1] <- "Numero_unique_DLoop"
colnames(d)[2] <- "Numero_unique_extrait"




# 2. Format columns contents ----------------------------------------------
# Everything allright with IDs columns: Numero_unique_D-Loop, Numero_unique_extrait, Numero_unique_specimen, Nom_Projet
# Also with Qualite_sequence, Sequence_consensus:Modifications

## 2.1. Responsable_Dloop -------------------------------------------------

table(d$Responsable_Dloop, useNA = "ifany")
d[d$Responsable_Dloop %in% "0.0", "Responsable_Dloop"] <- NA
d$Responsable_Dloop <- gsub("Sandrine et Frederique|Frederique et Sandrine|Sandrine, Frederique","Sandrine_Frederique",d$Responsable_Dloop)
d$Responsable_Dloop <- gsub("^Caroline$","CarolineChavarria",d$Responsable_Dloop)
d$Responsable_Dloop <- gsub(" |_","",d$Responsable_Dloop)
table(d$Responsable_Dloop, useNA = "ifany")


## 2.2. No_run ------------------------------------------------------------

table(d$No_run_F, useNA = 'ifany')
d[d$No_run_F %in% "0.0", "No_run_F"] <- NA  # NAs instead of 0.0
table(d$No_run_F, useNA = 'ifany')

table(d$No_run_R)
d[d$No_run_R %in% "0.0", "No_run_R"] <- NA  # NAs instead of 0.0
table(d$No_run_R, useNA = 'ifany')


## 2.3. No_plaque ---------------------------------------------------------

table(d$No_plaque_F, useNA = 'ifany')
d[d$No_plaque_F %in% "0.0", "No_plaque_F"] <- NA  # NAs instead of 0.0
d$No_plaque_F <- gsub("_$","",d$No_plaque_F)  # removes trailing underscore at the end of a string
d$No_plaque_F <- gsub("; |  et | et ","_",d$No_plaque_F)
d$No_plaque_F <- gsub(" ","",d$No_plaque_F)
table(d$No_plaque_F, useNA = 'ifany')

table(d$No_plaque_R)
d[d$No_plaque_R %in% "0.0", "No_plaque_R"] <- NA  # NAs instead of 0.0
d$No_plaque_R <- gsub("_$","",d$No_plaque_R)  # removes trailing underscore at the end of a string
d$No_plaque_R <- gsub("; |  et | et ","_",d$No_plaque_R)
d$No_plaque_R <- gsub(" ","",d$No_plaque_R)
table(d$No_plaque_R, useNA = 'ifany')


## 2.4. No_puits_F --------------------------------------------------------

table(d$No_puits_F, useNA = 'ifany')
d[d$No_puits_F %in% "0.0", "No_puits_F"] <- NA  # NAs instead of 0.0
table(d$No_puits_F, useNA = 'ifany')

table(d$No_puits_R, useNA = 'ifany')
d[d$No_puits_R %in% "0.0", "No_puits_R"] <- NA  # NAs instead of 0.0
table(d$No_puits_R, useNA = 'ifany')


## 2.5. Notes -------------------------------------------------------------

table(d$Notes, useNA = 'ifany')
d[d$Notes %in% "0.0", "Notes"] <- NA  # NAs instead of 0.0
d$Notes <- gsub("  "," ",d$Notes)
d$Notes <- gsub("[[:punct:][:blank:]]+", "_", d$Notes)  # https://stackoverflow.com/questions/29098801/removing-punctuations-from-text-using-r#comment46428660_29099172
d$Notes <- gsub("é", "e", d$Notes)  # there still are some accents here and there, remove them manually
d$Notes <- gsub("è", "e", d$Notes)  # there still are some accents here and there, remove them manually
d$Notes <- gsub("à", "a", d$Notes)  # there still are some accents here and there, remove them manually
d$Notes <- gsub("_$","",d$Notes)  # removes trailing underscore at the end of a string
table(d$Notes, useNA = 'ifany')




# 3. Save dataset ---------------------------------------------------------

write.csv(d, "Dloop_MOBELS.csv", row.names=F)  # Upload this directly on ACCESS file, D-Loop sheet


















