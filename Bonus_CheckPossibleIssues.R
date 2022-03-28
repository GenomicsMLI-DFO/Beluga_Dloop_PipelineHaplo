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


# Long haplotypes
table(d$N_nucl, useNA = "ifany")
table(h$N_nucl_615, useNA = "ifany")


which(h$Sequence_consensus %notin% d$Sequence_consensus)
which(d$Sequence_consensus %notin% h$Sequence_consensus)

table(is.na(d$Sequence_consensus))
table(is.na(h$Sequence_consensus))

table(d$haplotype, useNA = "ifany")
table(h$haplotype_615, useNA = "ifany")


# s <- read_excel("../ACCESS/20220318_MOBELS.xlsx", sheet = "Specimens", na = "NA")


























### 1.2.4. Format Sequence_consensus --------------------------------------

dt$Sequence_consensus <- toupper(dt$Sequence_consensus)  # Nucleotides in capital letters
dt$Sequence_consensus <- gsub("-", "", dt$Sequence_consensus)  # No breaks within sequences, remove '-' on the edges


#### 1.2.4.1. Filter by sequence length: remove short sequences -----------

dt$N_nucl2 <- as.integer(nchar(dt$Sequence_consensus))
table(dt$N_nucl2)
dt <- dt[!dt$N_nucl2 < 200,]  # remove any sequence that is very short (there should be none)

# Verify if and which ambiguities are included in sequences at present
nt <- paste(dt$Sequence_consensus, collapse = "")
table(strsplit(nt, split = ""))


#### 1.2.4.2. Identify duplicated sequences -------------------------------
# Necessary step to recognize different which sequences corresponds to which extraction/re-sequencing event in a fasta file

dloop <- dt %>%  # identify duplicated sequences by adding -2, -3, -4 after the ID of the specimen
  group_by(Numero_unique_specimen) %>%  # group by specimen ID
  mutate(Duplicated = rleid(Numero_unique_extrait)) %>%  # create new Duplicated column specifying which specimen is duplicated using numbers
  mutate(Numero_unique_specimen = paste(Numero_unique_specimen, Duplicated, sep = "-"))  # paste specimen ID with duplication number (add identified to Numero_unique_specimen)
dloop$Numero_unique_specimen <- gsub("-1", "", dloop$Numero_unique_specimen)  # unnecessary to specify which ones are unique or the first specimen of a series of duplicates









