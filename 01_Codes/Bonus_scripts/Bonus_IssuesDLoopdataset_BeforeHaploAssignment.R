# Info --------------------------------------------------------------------
#
# Overview: Issues in 'ACCESS' Beluga dataset (v. 20220603)
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2022-04-04
#


# Housekeeping ------------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
library(readxl)
library(dplyr)
# library(tidyr)

# Functions
"%nin%" <- Negate("%in%")




# Data --------------------------------------------------------------------


# Upload database

g <- data.frame(read_excel("../Documents/MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "Groupe", na = "NA"))
s <- data.frame(read_excel("../Documents/MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "Specimens", na = "NA"))
t <- data.frame(read_excel("../Documents/MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "Tissus", na = "NA"))
d <- data.frame(read_excel("../Documents/MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "D-Loop", na = "NA"))
sex <- data.frame(read_excel("../Documents/MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "05_qPCR", na = "NA"))

str(g)

str(s)

str(t)
# repeated column: Numero_unique_specimen
# 9 empty columns with no name
t <- t[,c(1:17)]
colnames(t)[2] <- "Numero_unique_specimen"

str(d)
# repeated column (three times): Numero_unique_extrait
# 36 empty columns before the last one
d <- d[,c(1:23)]
colnames(d)[2] <- "Numero_unique_extrait"

str(sex)
# 1 empty column
sex <- sex[,c(1:13,15:17)]


# Which specimens are in some sheets but not in other

## D-Loop not in Groupe
which(d$Numero_unique_specimen %nin% g$Numero_unique_reception_groupe)
d$Numero_unique_specimen[c(75,76,77)]
# "S_22_06671" "S_22_06673" "S_20_04348"

## D-Loop not in Specimens
which(d$Numero_unique_specimen %nin% s$Numero_unique_specimen)
d$Numero_unique_specimen[77]
# "S_20_04348"

## D-Loop not in Tissus
which(d$Numero_unique_specimen %nin% t$Numero_unique_specimen)
d$Numero_unique_specimen[77]
# "S_20_04348"


# Duplicated sequences
table(duplicated(d$Numero_unique_specimen))
dup <- d$Numero_unique_specimen[duplicated(d$Numero_unique_specimen)]
d1 <- d[d$Numero_unique_specimen %in% dup,]







