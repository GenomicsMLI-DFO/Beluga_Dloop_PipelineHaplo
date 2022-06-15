# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
# if(!require(tidyverse)){install.packages("tidyverse")}
library(readxl)
library(dplyr)

# if(!require(data.table)){install.packages("data.table")}
library(data.table)  # rleid function

# if(!require(stringr)){install.packages("stringr")}
library(stringr)  # str_count function

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings", force = TRUE)  # install Biostrings the first time you run this script
# BiocManager::install("msa", force = TRUE)  # install msa the first time you run this script
library(Biostrings)
library(msa)

# Functions
"%nin%" <- Negate("%in%")

library(stringr)


# 1. Data -----------------------------------------------------------------

## 1.1. Upload databases --------------------------------------------------

# Originally in ACCESS folder on Drive. Specify the path to the directory where the file is stored
d <- data.frame(read_excel("../../../Downloads/20220603_MOBELS.xlsx", sheet = "D-Loop", na = "NA"))  # remember to specify right path to beluga ACCESS dataset
# s <- read_excel("../../MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "Specimens", na = "NA")  # remember to specify right path to beluga ACCESS dataset
# g <- read_excel("../../MOBELS/DB/ACCESS/20220603_MOBELS.xlsx", sheet = "Groupe", na = "NA")  # remember to specify right path to beluga ACCESS dataset 

dnew <- data.frame(read_excel("../../MOBELS/DB/ACCESS/20220614_MOBELS_sans_doublons_DLoop_sexage.xlsx", sheet = "D-Loop", na = "NA"))

n1 <- read.csv2("../../../Downloads/20220330.csv", stringsAsFactors = F, )
n2 <- read.csv2("../../../Downloads/20220523.csv", stringsAsFactors = F)
n3 <- read.csv2("../../../Downloads/20220603_1.csv", stringsAsFactors = F)
n4 <- read.csv2("../../../Downloads/20220603_2.csv", stringsAsFactors = F)
n5 <- read.csv2("../../../Downloads/20220613.csv", stringsAsFactors = F)


## 1.2. Format input database for MSA -------------------------------------

### 1.2.1. Dloop ----------------------------------------------------------

str(d)  # 3993 rows
colnames(d)[2] <- "Numero_unique_extrait"

# Subset dataset: remove 'useless' columns
# d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
# d <- transform(d, Qualite_sequence = as.integer(d$Qualite_sequence),
#                N_nucl = as.integer(d$N_nucl))
d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, No_plaque_F, No_puits_F, No_plaque_R, No_puits_R, Sequence_consensus))

# Remove specimens without consensus sequence
d <- d[!is.na(d$Sequence_consensus),]  # removes 374 rows

str(dnew)  # 3993 rows
colnames(dnew)[2] <- "Numero_unique_extrait"

# Subset dataset: remove 'useless' columns
# d <- subset(d, select = c(Numero_unique_specimen, Numero_unique_extrait, Qualite_sequence, Sequence_consensus, N_nucl))
# d <- transform(d, Qualite_sequence = as.integer(d$Qualite_sequence),
#                N_nucl = as.integer(d$N_nucl))
dnew <- subset(dnew, select = c(Numero_unique_specimen, Numero_unique_extrait, No_plaque_F, No_puits_F, No_plaque_R, No_puits_R, Sequence_consensus))

# Remove specimens without consensus sequence
dnew <- dnew[!is.na(dnew$Sequence_consensus),]  # removes 374 rows



n <- rbind(n1,n2,n3,n4,n5)
table(n$ADN)
n$ADN <- gsub(".*\\ADN", "ADN", n$ADN)
n$ADN <- gsub("\\_F.*", "", n$ADN)
n$ADN <- gsub("\\_R.*", "", n$ADN)
table(nchar(n$ADN))
write.csv(n, "check.csv", row.names = F)

# ids <- unique(n$ADN)
# 
# d <- d[d$Numero_unique_extrait %in% ids,]
# dnew <- dnew[dnew$Numero_unique_extrait %in% ids,]

length(which(n$Seq %in% d$Sequence_consensus))
length(which(n$Seq %in% dnew$Sequence_consensus))

n$Length <- nchar(n$Seq)
n$Check_d <- ifelse(n$Seq %in% d$Sequence_consensus, T, F)
n$Check_dnew <- ifelse(n$Seq %in% dnew$Sequence_consensus, T, F)

for(i in 1:nrow(n)){
  dna <- n$ADN[i]
  n$Crsp <- ifelse(n$Seq[i] %in% d$Sequence_consensus[d$Numero_unique_extrait %in% dna], T, F)
}

for(i in 1:nrow(n)){
  dna <- n$ADN[i]
  n$Crsp_new <- ifelse(n$Seq[i] %in% dnew$Sequence_consensus[dnew$Numero_unique_extrait %in% dna], T, F)
}


dna <- "ADN_22_00225"

d$Sequence_consensus[d$Numero_unique_extrait %in% dna] == dnew$Sequence_consensus[dnew$Numero_unique_extrait %in% dna]




nn <- left_join(n, d[,c(2,7)], by = c("Seq"="Sequence_consensus"))
nn$Csp_old <- nn$ADN == nn$Numero_unique_extrait
nn <- arrange(nn, ADN)

nnn <- left_join(n, dnew[,c(2,7)], by = c("Seq"="Sequence_consensus"))
nnn$Csp_old <- nnn$ADN == nnn$Numero_unique_extrait
nnn <- arrange(nnn, ADN)

n_old <- n[n$Check_d %in% T, -5]
n_new <- n[n$Check_dnew %in% T, -4]

n_old1 <- left_join(n_old, d[,c(2,7)], by = c("Seq"="Sequence_consensus"))
n_old1$Crsp <- n_old1$ADN == n_old1$Numero_unique_extrait
table(n_old1$Crsp)
n_old1 <- arrange(n_old1, ADN)

n$Check_d <- ifelse(n$Seq %in% d$Sequence_consensus, d$Numero_unique_extrait[d$Sequence_consensus %in% n$Seq], "F")
n$Check_dnew <- ifelse(n$Seq %in% dnew$Sequence_consensus, dnew$Numero_unique_extrait[dnew$Sequence_consensus %in% n$Seq], "F")
n$ids_d <- n$ADN == n$Check_d
n$ids_dnew <- n$ADN == n$Check_dnew


which(d$Sequence_consensus %in% n$Seq[331])


















