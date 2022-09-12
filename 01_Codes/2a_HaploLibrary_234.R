# Info --------------------------------------------------------------------
# 
# Authors: Benjamin Hornoy, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-12-17
# 
# Overview: Compiling the D-loop haplotype library
# 
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
if(!require(stringr)){install.packages("stringr")}
library(stringr)

library(adegenet)

# Functions





# 1. Data -----------------------------------------------------------------

data <- read.table("00_Data/02_dloop_clean/Sequences_Dloop234_n3441.txt", header = T)
str(data)
colnames(data)[2] <- "seq"




# 2. Sequences' quality ---------------------------------------------------

sequences <- toupper(data$seq)
nucleotides <- c("A","T","C","G")
ambig_nucl <- c("N","R","Y","K","M","S","W","B","D","H","V")
seq_len_exp <- 234  # expected sequence length

res <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(res) <- c("N.nucl","N.ATCG", "N.ambig", "N.manquants")
for (i in 1:length(sequences)){
    if(is.na(sequences[i])){
        res[i,] <- c("NA","NA","NA","NA")	
    }else{
        seq_len <- sum(str_count(sequences[i],c(nucleotides, ambig_nucl)))
    }	
    res[i,] <- c(seq_len, sum(str_count(sequences[i], nucleotides)), sum(str_count(sequences[i], ambig_nucl)), seq_len_exp-seq_len)
}
data1 <- cbind(data, res)




# 3. Sequence: usable or not? ---------------------------------------------
# Considered usable if 100% complete on the minimal sequence AND no ambiguous nucleotides in within the minimal sequence

dna <- fasta2DNAbin("00_Data/01_fasta/Beluga_234bp_n3441.fasta")  # import fasta sequence
dna <- DNAbin2genind(dna, polyThres=0)  # trasform DNAbin object into genind object
dna  # for info
snpPos <- locNames(dna)  # vector with position of polymorphisms (SNPs) within sequences

# obtain info on minimal sequence
seq_start <- as.numeric(snpPos[1])  # position of first polymorphic site (6)
seq_stop <- as.numeric(tail(snpPos, 1))  # position of last polymorphic site (230)

util <- data.frame(matrix(ncol=1, nrow=0))
colnames(util) <- c("seq_utilisable")

for (i in 1:length(sequences)){
    if (is.na(sequences[i])){
        util[i,] <- 0
    } else if (sum(str_count(substr(sequences[i],seq_start,seq_stop),nucleotides)) == seq_stop-seq_start+1) {
        util[i,] <- 1  # sum no of ACTG should = minimal sequence (225 nt). Rationale: if ambiguities are present sum of A,C,T,G < 570
    }else{
        util[i,] <- 0
    }
}
data2 <- cbind(data1, util)
write.csv(data2, "00_Data/02_dloop_clean/Sequences_Dloop234_n3441.csv", row.names = F)




# # 4. Compile haplotype library --------------------------------------------
# # Use this section to create 'original' haplotype library
# # If haplotype library is already present, use 5. Extend haplotype library below
# 
# lib <- data.frame(matrix(nrow = 0, ncol = 2))
# colnames(lib) <- c("hapl","seq")
# lib[1,] <- c("HS001", data2$seq[1])
# 
# nlib <- data.frame(matrix(ncol=2, nrow=0))
# a <- length(lib$seq)  # nb haplo in starting library
# for(i in 1:length(data2$seq)){
#     if(data2$seq_utilisable[i] == 1) {
#         if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {
# 
#         } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {
# 
#         } else {
#             a=a+1
#             if(a < 10) {
#                 nlib <- rbind(nlib, cbind(paste("HS00",a,sep=""), data2$seq[i]))
#             } else if(a < 100 & a >= 10) {
#                 nlib <- rbind(nlib, cbind(paste("HS0",a,sep=""), data2$seq[i]))
#             } else {
#                 nlib <- rbind(nlib,cbind(paste("HS",a,sep=""), data2$seq[i]))
#             }
#         }
#     } else {
#     }
# }
# colnames(nlib) <- c("hapl","seq")
# lib_fin <- rbind(lib,nlib)
# 
# 
# ## 4.1. Save haplotype library --------------------------------------------
# 
# ### 4.1.1. Create directory -----------------------------------------------
# 
# if(dir.exists("./02_Results/00_libraries") == F){
#     dir.create("./02_Results/00_libraries")
# }
# 
# ### 4.1.2. Save haplotype library ----------------------------------------
# 
# write.csv(lib_fin, file = paste("02_Results/00_libraries/","librairie_", length(lib_fin$hapl), "_haplotypes234.csv", sep=""), row.names = F)




# 5. Extend haplotype library ---------------------------------------------

lib <- read.csv("02_Results/00_libraries/librairie_53_haplotypes234.csv")  # upload short haplo library
colnames(lib) <- c("hapl","seq")  # if it's not already the case
table(nchar(lib$seq))  # all haplotypes are 234 nt long


## 5.1. Detect new haplotypes --------------------------------------------

a <- length(lib$seq) # nb haplo in starting library
nlib <- data.frame(matrix(ncol=2, nrow=0))  # creates a new library

for(i in 1:length(data2$seq)){  # compile new library
    if(data2$seq_utilisable[i] == 1) {
        if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {

        } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {

        } else {
            a=a+1
            if(a < 100) {
                nlib <- rbind(nlib, cbind(paste("HS0",a,sep=""), data2$seq[i]))
            } else {
                nlib <- rbind(nlib,cbind(paste("HS",a,sep=""), data2$seq[i]))
            }
        }
    } else {
    }
}  # library will be empty if no new haplotypes are found
colnames(nlib) <- c("hapl","seq")
lib_fin <- rbind(lib,nlib)

# Save the updated haplotype library
write.csv(lib_fin, file = paste("02_Results/00_libraries/","librairie_", length(lib_fin$hapl), "_haplotypes234.csv", sep=""), row.names = F)


### WARNING! Sequences might contain ambiguous or empty sites (out of minimal sequence), which could produce problems in the future once these sites will be detected
### as polymorphisms. We must either authomatize or check by hand that new sequences introduced in the library won't contain ambiguities or empty sites. Hopefully, no
### new haplotypes will be represented by one individual (one sequence) with ambiguities or empty sites...


## 5.2. New library: ambiguities? ----------------------------------------

ambigSites <- c(ambig_nucl, "-")

for (i in 1:length(ambigSites)){
    if (length(grep(ambigSites[i], lib_fin$seq))==0){
    }else{
        print(paste(ambigSites[i]," a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL",grep(ambigSites[i], lib_fin$seq),sep=""))
    }
}
# if no text is printed, everything is fine
