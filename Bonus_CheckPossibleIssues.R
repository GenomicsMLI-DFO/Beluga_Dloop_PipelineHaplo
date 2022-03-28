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


str(d)  # 3643 rows
str(h)  # 3643 rows

h <- arrange(h, Numero_unique_extrait)
d <- arrange(d, Numero_unique_extrait)

# Subset dataset: remove 'useless' columns
d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl, haplotype))
h <- subset(h, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl_234, haplotype_234, N_nucl_615, haplotype_615))

d <- transform(d, Qualite_sequence = as.integer(d$Qualite_sequence),
               N_nucl = as.integer(d$N_nucl))
h <- transform(h, Qualite_sequence = as.integer(d$Qualite_sequence),
               N_nucl_234 = as.integer(h$N_nucl_234),
               N_nucl_615 = as.integer(h$N_nucl_615))

# Remove specimens without consensus sequence
table(is.na(d$Sequence_consensus))
table(is.na(h$Sequence_consensus))
d <- d[!is.na(d$Sequence_consensus),]  # removes 315 rows
h <- h[!is.na(h$Sequence_consensus),]  # removes 315 rows


# Take a look at sequence quality - defined by lab technicians (?)
table(d$Qualite_sequence, useNA = 'ifany')
#   1    2    3    4   11 
#3300   17    5    1    5
table(h$Qualite_sequence, useNA = 'ifany')
#    1    2    3    4   11 
# 3300   17    5    1    5


# Number of nucleotides
table(d$N_nucl, useNA = "ifany")
table(h$N_nucl_615, useNA = "ifany")

# How many long haplo are assigned
table(d$haplotype, useNA = "ifany")
table(h$haplotype_615, useNA = "ifany")


table(h$Numero_unique_extrait == d$Numero_unique_extrait, useNA = "ifany")
table(h$Sequence_consensus == d$Sequence_consensus, useNA = "ifany")
table(h$haplotype_615 == d$haplotype, useNA = "ifany")

dloop <- left_join(d[,names(d) %in% c("Numero_unique_specimen","Numero_unique_extrait","Sequence_consensus")],
                   h[,names(h) %in% c("Numero_unique_specimen","Numero_unique_extrait","Sequence_consensus")],
                   by = c("Numero_unique_specimen", "Numero_unique_extrait"))
dloop$same <- dloop$Sequence_consensus.x == dloop$Sequence_consensus.y
table(dloop$same)


