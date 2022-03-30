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
library(readxl)
library(dplyr)
# library(data.table)  # rleid function

# Functions
"%notin%" <- Negate("%in%")



# 1. Data -----------------------------------------------------------------

## 1.1. Upload databases --------------------------------------------------

# Originally in ACCESS folder on Drive. Specify the path to the directory where the file is stored
h <- read.csv("Dloop_haplo_n3643.csv")
d <- read.csv("Dloop_MOBELS.csv")
s <- read_excel("../ACCESS/20220328_MOBELS.xlsx", sheet = "Specimens", na = "NA")

str(d)  # 3643 rows
str(h)  # 3643 rows
str(s)

h <- arrange(h, Numero_unique_extrait)
d <- arrange(d, Numero_unique_extrait)
s <- arrange(s, Numero_unique_specimen)
colnames(s)[29] <- "Extended_Type"

# Subset dataset: remove 'useless' columns
d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl, haplotype))
h <- subset(h, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl_615, haplotype_615))
s <- subset(s, select = c(Numero_unique_specimen, Nom_commun, Short_Haplo, Extended_Type))
table(is.na(s$Extended_Type))
table(is.na(s$Short_Haplo))

d <- transform(d, Qualite_sequence = as.integer(d$Qualite_sequence),
               N_nucl = as.integer(d$N_nucl))
h <- transform(h, Qualite_sequence = as.integer(d$Qualite_sequence),
               N_nucl_615 = as.integer(h$N_nucl_615))

dt <- left_join(d, h[,names(h) %notin% c("Qualite_sequence", "Sequence_consensus")], by = c("Numero_unique_specimen","Numero_unique_extrait"))
# dt$Same <- dt$Sequence_consensus.x == dt$Sequence_consensus.y
# table(dt$Same, useNA = "ifany")  # T = 3328 F = 315, all good
dt <- left_join(dt, s, by = "Numero_unique_specimen")

# Haplotypes
table(dt$haplotype, useNA = "ifany")  # there are three 0 that should be changed to NA
table(dt$haplotype_615, useNA = "ifany")
dt$haplotype[dt$haplotype %in% 0] <- NA


# Compare which specimens did not have an haplo assigned
table(is.na(dt$haplotype))  # F = 3296; T = 347
table(is.na(dt$haplotype_615))  # F = 3273; T = 370

# less specimes with unassigned haplo in 'original' D-Loop excel sheet: why?
nonew <- which(dt$Numero_unique_specimen[is.na(dt$haplotype_615)] %notin% dt$Numero_unique_specimen[is.na(dt$haplotype)])  # 28 specimens in h1 do not have an haplo but have it in d1
# 7   8  14  15  17  18  20  23  30  34  40  86 143 167 173 255 258 269 270 271 272 274 303 308 353 354 356 357 OF THE IS.NA LIST (not comprehensive dt)
noold <- which(dt$Numero_unique_specimen[is.na(dt$haplotype)] %notin% dt$Numero_unique_specimen[is.na(dt$haplotype_615)])  # 7 specimens in d1 do not have an haplo but have it in h1
# 9  92 106 234 283 284 328 OF THE IS.NA LIST (not comprehensive dt)
nonew_id <- dt$Numero_unique_specimen[is.na(dt$haplotype_615)]
nonew_id <- nonew_id[nonew]
noold_id <- dt$Numero_unique_specimen[is.na(dt$haplotype)]
noold_id <- noold_id[noold]
dt$MissingNew <- ifelse(dt$Numero_unique_specimen %in% nonew_id, 1, 0)
dt$MissingOld <- ifelse(dt$Numero_unique_specimen %in% noold_id, 1, 0)


# Despite some sequences are seemingly long enough (nt > 59 - minimal sequence length defined in 1_SNPsMinimalSequence.R, e.g. S_20_1034), 597 nt are necessary WITHIN
# the most extreme SNPs found (nt positions: 15 and 611 as of 20220329). A few sequences ends too early despite starting with the classic ACTACG sequence.
# The minimal sequence length and extreme boundaries have changed since Frederique's last haplo assignment (more specimens haplotyped)
# You cab check this using 2b_HaploLibrary_615.R, LL 43 - 90

# Verify if sequence length is the same of old dataset (after removing "---")
d$Sequence_consensus <- gsub("-", "", d$Sequence_consensus)
sequences <- left_join(subset(d, select = -c(Qualite_sequence, N_nucl, haplotype)),
                       subset(h, select = -c(Qualite_sequence, N_nucl_615, haplotype_615)),
                       by = c("Numero_unique_specimen","Numero_unique_extrait"))
table(sequences$Sequence_consensus.x == sequences$Sequence_consensus.y)  # ALL GOOD

x <- read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "D-Loop", na = "NA")  # new excel file, to see if everything's fine there as well - change name every time the
# D-Loop file is updated
colnames(x)[2] <- "Numero_unique_extrait"
x <- subset(x, select = c(Numero_unique_specimen,Numero_unique_extrait,Sequence_consensus))
sequences <- left_join(sequences, x, by = c("Numero_unique_specimen","Numero_unique_extrait"))
table(sequences$Sequence_consensus.x == sequences$Sequence_consensus)  # ALL GOOD
table(sequences$Sequence_consensus.y == sequences$Sequence_consensus)  # ALL GOOD
