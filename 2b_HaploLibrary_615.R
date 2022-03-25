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
# if(!require(stringr)){install.packages("stringr")}
library(stringr)

library(adegenet)

# Functions





# 1. Data -----------------------------------------------------------------

data <- read.table("Sequences_Dloop615_n3314.txt", header = T)
str(data)
colnames(data)[2] <- "seq"




# 2. Sequences' quality ---------------------------------------------------

sequences <- toupper(data$seq)
nucleotides <- c("A","T","C","G")
ambig_nucl <- c("N","R","Y","K","M","S","W","B","D","H","V")
seq_len_exp <- 615  # expected sequence length

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

dna <- fasta2DNAbin("fasta/Beluga_615bp_n3314.fasta")  # import fasta sequence
dna <- DNAbin2genind(dna, polyThres=0)  # trasform DNAbin object into genind object
dna  # for info
snpPos <- locNames(dna)  # vector with position of polimorphisms (SNPs) within sequences

# obtain info on minimal sequence
seq_start <- as.numeric(snpPos[1])  # position of first polymorphic site (15)
seq_stop <- as.numeric(tail(snpPos, 1))  # position of last polymorphic site (611)

util <- data.frame(matrix(ncol=1, nrow=0))
colnames(util) <- c("seq_utilisable")

for (i in 1:length(sequences)){
  if (is.na(sequences[i])){
    util[i,] <- "no"
  } else if (sum(str_count(substr(sequences[i],seq_start,seq_stop),nucleotides)) == seq_stop-seq_start+1) {
    util[i,] <- "yes"  # sum no of ACTG should = minimal sequence (225 nt). Rationale: if ambiguities are present sum of A,C,T,G < 570
  }else{
    util[i,] <- "no"
  }
}
data2 <- cbind(data1, util)
write.csv(data2, "Sequences_Dloop615_n3314.csv", row.names = F)




# # 4. Compile haplotype library --------------------------------------------
# # Use this section to create 'original' haplotype library
# # If haplotype library is already present, use 5. Extend haplotype library below
# 
# lib <- data.frame(matrix(nrow = 0, ncol = 2))
# colnames(lib) <- c("hapl","seq")
# lib[1,] <- c("HL001", data2$seq[1])
# 
# nlib <- data.frame(matrix(ncol=2, nrow=0))
# a <- length(lib$seq)  # nb haplo in starting library
# for(i in 1:length(data2$seq)){
#     if(data2$seq_utilisable[i] == "yes") {
#         if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {
# 
#         } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {
# 
#         } else {
#             a=a+1
#             if(a < 10) {
#                 nlib <- rbind(nlib, cbind(paste("HL00",a,sep=""), data2$seq[i]))
#             } else if(a < 100 & a >= 10) {
#                 nlib <- rbind(nlib, cbind(paste("HL0",a,sep=""), data2$seq[i]))
#             } else {
#                 nlib <- rbind(nlib,cbind(paste("HL",a,sep=""), data2$seq[i]))
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
# projdir <- getwd()
# outdir <- file.path(projdir, "libraries")
# dir.create(outdir)
# 
# ### 4.1.2. Save haplotype library ----------------------------------------
# 
# write.csv(lib_fin, file = paste("libraries/","librairie_", length(lib_fin$hapl), "_haplotypes615.csv", sep=""), row.names = F)




# 5. Extend haplotype library ---------------------------------------------

lib <- read.csv("libraries/librairie_137_haplotypes615.csv")  # upload short haplo library
colnames(lib) <- c("hapl","seq")  # if it's not already the case
table(nchar(lib$seq))  # all haplotypes are 615 nt long


## 5.1. Detect new haplotypes --------------------------------------------

a <- length(lib$seq) # nb haplo in starting library
nlib <- data.frame(matrix(ncol=2, nrow=0))  # creates a new library

for(i in 1:length(data2$seq)){  # compile new library
  if(data2$seq_utilisable[i] == "yes") {
    if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {
      
    } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {
      
    } else {
      a=a+1
      if(a < 100) {
        nlib <- rbind(nlib, cbind(paste("HL0",a,sep=""), data2$seq[i]))
      } else {
        nlib <- rbind(nlib,cbind(paste("HL",a,sep=""), data2$seq[i]))
      }
    }
  } else {
  }
}  # library will be empty if no new haplotypes are found
colnames(nlib) <- c("hapl","seq")
lib_fin <- rbind(lib,nlib)

# Save the updated haplotype library
write.csv(lib_fin, file = paste("libraries/","librairie_", length(lib_fin$hapl), "_haplotypes615.csv", sep=""), row.names = F)


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
# [1] "R a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL126"
# [1] "- a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL104"

# Use minimal sequence to understand if multiple sequences (use minimal sequence) have the same haplotype
# data2[substr(data2$seq, seq_start, seq_stop) %in% substr(lib_fin[lib_fin$hapl %in% "HL126","seq"], seq_start, seq_stop),]
# HL126 found in S_20_03408, S_20_03546. S_20_03546 with A instead of R
data2[substr(data2$seq, seq_start, seq_stop) %in% substr(lib_fin[lib_fin$hapl %in% "HL104","seq"], seq_start, seq_stop),] 
# HL104 only found for S_20_01854
# data2[substr(data2$seq, seq_start, seq_stop)==substr(lib_fin[lib_fin$hapl=="HL068","seq"], seq_start, seq_stop),]
# HL068 also found in S_20_01102-2, S_20_01436, S_20_02110, S_20_02167, S_20_02568, S_20_03122, S_20_03198, S_20_03652, S_20_03679 and with C instead of S


## 5.3. Update haplo for sequences without ambiguities -------------------
# Change sequences in lib haplo with ambiguities (outside minimal sequence) with 'clean' sequence

# lib_fin[lib_fin$hapl=="HL068", "seq"] <- data2$seq[data2$ID=="S_20_01102-2"]
# lib_fin[lib_fin$hapl %in% "HL126", "seq"] <- data2$seq[data2$ID %in% "S_20_03546"]

## 5.4. Write final library ----------------------------------------------

write.csv(lib_fin, file = paste("libraries/","librairie_", length(lib_fin$hapl), "_haplotypes615.csv", sep=""), row.names = F)




