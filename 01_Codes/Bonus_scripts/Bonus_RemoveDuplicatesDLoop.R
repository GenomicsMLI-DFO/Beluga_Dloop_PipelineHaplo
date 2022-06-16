# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-12-17
# 
# Overview: D-Loop sheet has duplicated specimens (all with NA haplo - keep one; some with NA and some with HL; all with HL; etc)
# Removes 'useless' duplicated (haplo NAs when HL is known, duplicated HL)
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

d <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220616_MOBELS.xlsx", sheet = "D-Loop", na = "NA"))
s <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220616_MOBELS.xlsx", sheet = "Specimens", na = "NA"))

str(d)
colnames(d)[c(2,24)] <- c("Numero_unique_extrait","Numero_unique_extrait_2")
table(duplicated(d$Numero_unique_specimen))  # 3683 rows expected after removing all duplicates using script below
# FALSE  TRUE 
#  3683   244

str(s)
s <- s[, c("Numero_unique_specimen", "Nom_commun")]

d <- left_join(d, s, by = "Numero_unique_specimen")
table(d$Nom_commun, useNA = "ifany")
# Beluga Beluga_Narval  Narval 
#   3909             3      15



# 2. Find duplicated specimens --------------------------------------------

nbel <- d[d$Nom_commun %nin% "Beluga",] %>% 
  subset(select = -Nom_commun)
bel <- d[d$Nom_commun %in% "Beluga",]

table(duplicated(bel$Numero_unique_specimen))  # 243 duplicated specimens
id <- bel$Numero_unique_specimen[duplicated(bel$Numero_unique_specimen)]

## 2.1. Create duplicated and unique specimen datasets --------------------

uni <- bel[bel$Numero_unique_specimen %nin% id, ] %>% 
  subset(select = -Nom_commun)
table(duplicated(uni$Numero_unique_specimen))
dup <- bel[bel$Numero_unique_specimen %in% id, ]
length(unique(dup$Numero_unique_specimen))  # 205 - this should also be the number of rows once I removed duplicates
# dup <- dup %>%
#   group_by(Numero_unique_specimen) %>% 
#   mutate(Dup = length(Numero_unique_specimen))
# table(dup$Dup)  # to confirm that specimens here are all duplicated
# #   2   3   4   5 
# # 344  87  12   5



# 3. Remove duplicates ----------------------------------------------------

## 3.1. Identify which duplicates have or don't have HL (and HS) ----------

dup <- dup %>% 
  arrange(Numero_unique_specimen, Numero_unique_extrait) %>% 
  mutate(Nchar = nchar(Sequence_consensus)) %>% 
  group_by(Numero_unique_specimen) %>% 
  mutate(keep = ifelse(all(is.na(haplotype_615)), 1, 0)) %>%  # identify if all duplicates of a specimen have no HL (1) or if some have been haplotyped (0)
  mutate(keep234 = ifelse(keep %in% 1 & all(is.na(haplotype_234)), 1, 0)) %>%  # identify if duplicates of a specimen without HL have (0) or not (1) the HS
  data.frame()

dup1 <- dup[dup$keep %in% 0, ]  # duplicates with HL
dup2 <- dup[dup$keep %in% 1 & dup$keep234 %in% 0, ]  # duplicates without HL but with HS
dup3 <- dup[dup$keep %in% 1 & dup$keep234 %in% 1, ]  # duplicates without HL and HS

# Verify if ID in dup1 are not in dup2 or dup3, and dup2 are not in dup3
which(dup1$Numero_unique_specimen %in% dup2$Numero_unique_specimen)
which(dup1$Numero_unique_specimen %in% dup3$Numero_unique_specimen)
which(dup2$Numero_unique_specimen %in% dup3$Numero_unique_specimen)


## 3.2. Duplicated with HL ------------------------------------------------

### 3.2.1. Remove specimens without HL ------------------------------------

table(is.na(dup1$haplotype_615))
dup1 <- dup1[!is.na(dup1$haplotype_615), ]  # remove specimens without haplo

### 3.2.2. Identify individuals that are still duplicated -----------------

# dup1 <- dup1 %>%
#   group_by(Numero_unique_specimen) %>%
#   mutate(Dup = length(Numero_unique_specimen)) %>% 
#   data.frame()
# table(dup1$Dup)  # to confirm if some specimens are now unique
# #  1   2   3   4 
# # 70 208  51   4

ids <- dup1$Numero_unique_specimen[duplicated(dup1$Numero_unique_specimen)]
dup1.1 <- dup1[dup1$Numero_unique_specimen %in% ids,]  # still duplicated
dup1.2 <- dup1[dup1$Numero_unique_specimen %nin% ids,]  # not duplicated anymore

### 3.2.3. Removes seq < 615 nt or with ambiguities -----------------------

dup1.1 <- dup1.1[!(dup1.1$N_ambig_615 > 0 | dup1.1$N_manquants_615 > 0),]

### 3.2.4. Identify individuals that are still duplicated -----------------

ids <- dup1.1$Numero_unique_specimen[duplicated(dup1.1$Numero_unique_specimen)]
dup1.1.1 <- dup1.1[dup1.1$Numero_unique_specimen %in% ids,]  # still duplicated
dup1.1.2 <- dup1.1[dup1.1$Numero_unique_specimen %nin% ids,]  # not duplicated anymore

### 3.2.5. Duplicates: remove those with shortest consensus sequence ------

dup1.1.1 <- dup1.1.1 %>% 
  arrange(Numero_unique_specimen, Nchar)
dup1.1.1 <- dup1.1.1[!duplicated(dup1.1.1$Numero_unique_specimen, fromLast = T),]

dup1 <- rbind(dup1.1.1, dup1.1.2, dup1.2)
col_to_remove <- c("Nom_commun","Nchar","keep","keep234")
dup1 <- dup1[, !names(dup1) %in% col_to_remove]

## 3.3. Duplicated without HL but with HS ---------------------------------

dup2 <- dup2[!is.na(dup2$haplotype_234), ]
dup2 <- dup2[, !names(dup2) %in% col_to_remove]


## 3.4. Duplicated without HL nor HS --------------------------------------

dup3 <- dup3 %>%
  group_by(Numero_unique_specimen) %>% 
  mutate(keep1 = ifelse(all(is.na(Nchar)), 1, 0))

dup3.1 <- dup3[dup3$keep1 %in% 1, ]  # duplicates without consensus sequence
dup3.2 <- dup3[dup3$keep1 %in% 0, ]  # duplicates with consensus sequence

col_to_remove <- c("Nom_commun","Nchar","keep","keep234","keep1")

dup3.1 <- dup3.1[!duplicated(dup3.1$Numero_unique_specimen), !names(dup3.1) %in% col_to_remove]  # doesn't matter which one is removes
dup3.2 <- dup3.2[!is.na(dup3.2$Nchar), !names(dup3.2) %in% col_to_remove]

dup3 <- rbind(dup3.1, dup3.2)

dup <- rbind(dup1, dup2, dup3)  # 205 rows, as it should given l64



# 4. Clean D-Loop file ----------------------------------------------------

dt <- rbind(dup, uni, nbel)
write.csv(dt, "./02_Results/02_ACCESS/Dloop_haplo_n3684_no_duplicates.csv")


table(duplicated(dt$Numero_unique_specimen))

table(duplicated(d$Numero_unique_specimen))


