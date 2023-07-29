# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2022-06-16
# 
# Overview:  of specimens that should be resequenced to get long haplotype
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
"%nin%" <- Negate("%in%")



# 1. Data -----------------------------------------------------------------

d <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220617_MOBELS.xlsx", sheet = "D-Loop", na = "NA"))
s <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220617_MOBELS.xlsx", sheet = "Specimens", na = "NA"))
g <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220617_MOBELS.xlsx", sheet = "Groupe", na = "NA"))

str(s)
s <- s[, c("Numero_unique_specimen","Nom_commun")]

str(d)
d <- subset(d, select = c(Numero_unique_specimen,Numero_unique_extrait,Sequence_consensus,N_nucl_615,N_ambig_615,N_manquants_615,haplotype_615))

str(g)
g <- subset(g, select = c(Numero_unique_reception_groupe,Region_echantillonnage,Lieu_echantillonnage,Annee_echantillonnage:Jour_echantillonnage))

dt <- left_join(d, s, by = "Numero_unique_specimen")
dt <- left_join(dt, g, by = c("Numero_unique_specimen"="Numero_unique_reception_groupe"))




# 2. To be resampled from 2021 hunt ---------------------------------------

table(dt$Annee_echantillonnage, useNA = "ifany")
nrow(dt[grepl("S_22_", dt$Numero_unique_specimen),])

dt <- dt[dt$Annee_echantillonnage %in% 2021 | grepl("S_22_", dt$Numero_unique_specimen),]
table(dt$haplotype_615, useNA = "ifany")

dt$Numero_unique_specimen[is.na(dt$haplotype_615)]











dt <- dt[!is.na(dt$Region_echantillonnage),]
dt <- dt[!is.na(dt$Mois_echantillonnage),]
dt <- dt[dt$Nom_commun %in% "Beluga",]  # keep beluga only

# which specimens actually have an haplo because they have been sequences multiple times
rem <- which(dt$Numero_unique_specimen %in% d$Numero_unique_specimen[!is.na(d$haplotype_615)])
dt <- dt[-rem,]  # 124 specimens that are extracted once (or even multiple times) but do not have any haplotype associated with them

# nchar(d[d$Numero_unique_specimen %in% "S_20_03210", "Sequence_consensus"])

# # Prepare list of samples to resamples with info on where to find them
# ext <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Extraits_ADN_ARN", na = "NA"))
# tis <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Tissus", na = "NA"))
# 
# str(ext)
# colnames(ext)[2:4] <- c("Numero_unique_specimen","Numero_unique_tissu","Numero_unique_extrait")
# ext <- subset(ext, select = c(Numero_unique_specimen,Numero_unique_tissu,Numero_unique_extrait,Numero_boite_extrait))
# 
# str(tis)
# colnames(tis)[2] <- "Numero_unique_specimen"
# tis <- subset(tis, select = c(Numero_unique_specimen,Numero_unique_tissu,No_Local_entreposage_tissu,No_congelateur_frigo,No_boite_entreposage_tissu))
# 
# samples <- left_join(dt, tis, by = c("Numero_unique_specimen"))






