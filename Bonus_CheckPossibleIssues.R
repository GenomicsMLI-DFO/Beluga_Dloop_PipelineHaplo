# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2022-03-28
# 
# Overview: Verify if you introduce any possible mistake in the D-Loop dataset after running the pipeline
# Qualitative check
# 
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
# library(readxl)
library(dplyr)
# library(data.table)  # rleid function

# Functions
"%notin%" <- Negate("%in%")



# 1. Data -----------------------------------------------------------------

## 1.1. Upload databases --------------------------------------------------

# Originally in ACCESS folder on Drive. Specify the path to the directory where the file is stored
h <- read.csv("Dloop_haplo_n3643.csv")
d <- read.csv("Dloop_MOBELS.csv")

str(h)

h <- arrange(h, Numero_unique_extrait)
d <- arrange(d, Numero_unique_extrait)

which(h$Sequence_consensus %notin% d$Sequence_consensus)
which(d$Sequence_consensus %notin% h$Sequence_consensus)

table(is.na(d$Sequence_consensus))
table(is.na(h$Sequence_consensus))

table(d$haplotype, useNA = "ifany")
table(h$haplotype_615, useNA = "ifany")


# s <- read_excel("../ACCESS/20220318_MOBELS.xlsx", sheet = "Specimens", na = "NA")


























