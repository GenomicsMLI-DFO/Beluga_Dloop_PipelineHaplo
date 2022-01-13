# Info --------------------------------------------------------------------

# Compiling the d_haplotypes library
# Short (234 bp) haplotypes
# 
# Benjamin Hornoy
# Minor changes by Luca Montana
#
#

getwd()
rm(list = ls())

# Libraries and function --------------------------------------------------
if(!require(stringr)){install.packages("stringr")}
library(stringr)

library(adegenet)


# Data --------------------------------------------------------------------
data <- read.table("Sequences_Dloop234_n3314.txt", header = T)
str(data)
colnames(data)[2] <- "seq"

# Include info on quality of sequences ------------------------------------
sequences <- toupper(data$seq)
nucleotides <- c("A","T","C","G")
ambig_nucl <- c("N","R","Y","K","M","S","W","B","D","H","V")
seq_len_exp <- 234  # valeur de la longueur attendue de la séquence

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


# Include columns specifying if sequence if usable ------------------------
# Considered usable if 100% complete on the minimal sequence
# AND no ambiguous nucleotides in within the minimal sequence
dna <- fasta2DNAbin("Beluga_234bp_n3314.fasta")
dna <- DNAbin2genind(dna, polyThres=0)
dna  # for info
snpPos <- locNames(dna)

# obtain info on minimal sequence
seq_start <- as.numeric(snpPos[1])  # position of first polymorphic site (6)
seq_stop <- as.numeric(tail(snpPos, 1))  # position of last polymorphic site (234)

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
write.csv(data2, "Sequences_Dloop234_n3314.csv", row.names = F)


# Compile new haplotype library -------------------------------------------
lib <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(lib) <- c("hapl","seq")
lib[1,] <- c("HS001", data2$seq[1])

nlib <- data.frame(matrix(ncol=2, nrow=0))
a <- length(lib$seq)  # nb haplo in starting library
for(i in 1:length(data2$seq)){
    if(data2$seq_utilisable[i] == "yes") {
        if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {

        } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {

        } else {
            a=a+1
            if(a < 10) {
                nlib <- rbind(nlib, cbind(paste("HS00",a,sep=""), data2$seq[i]))
            } else if(a < 100 & a >= 10) {
                nlib <- rbind(nlib, cbind(paste("HS0",a,sep=""), data2$seq[i]))
            } else {
                nlib <- rbind(nlib,cbind(paste("HS",a,sep=""), data2$seq[i]))
            }
        }
    } else {
    }
}
colnames(nlib) <- c("hapl","seq")
lib_fin <- rbind(lib,nlib)

# Save the haplotype library
write.csv(lib_fin, file=paste("librairie_", length(lib_fin$hapl), "_haplotypes234.csv", sep=""), row.names=F)
# Code above used to create new 'official' HS library
# Next time new sequences are added, use code below : Expand haplotype library


# Expand haplotype library ------------------------------------------------
# lib <- read.csv("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/dloop/librairie_139_haplotypes.csv", )
# colnames(lib) <- c("hapl","seq") # if it's not already the case
# table(nchar(lib$seq))  # 255 = 1; 615 = 138: there is one haplo identified on 255 nt
# 
# # detect new haplotypes
# a <- length(lib$seq) # nb haplo in starting library
# nlib <- data.frame(matrix(ncol=2, nrow=0))
# 
# for(i in 1:length(data2$seq)){
#     if(data2$seq_utilisable[i] == "yes") {
#         if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)) {
# 
#         } else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)) {
# 
#         } else {
#             a=a+1
#             if(a < 100) {
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

# # Save the updated haplotype library
# write.csv(lib_fin, file=paste("librairie_", length(lib_fin$hapl), "_haplotypes615.csv", sep=""), row.names=F)


### ATTENTION! La "seq" peut contenir des sites ambigus ou manquants (hors de la séquence minimale), ce qui pourra poser pb dans le futur si ces sites sont détectés comme polymorphes dans le futur.
### Il faut donc automatiser ou s'assurer à la mitaine que la séquence mise dans la librairie ne contient pas d'ambigus, en espérant qu'aucun haplotype ne soit représenté que par un seul individu, qui a des ambigus...

# vérifier si des ambigus sont présents dans la nouvelle librairie (si oui, voir ci-dessus...)
ambigSites <- c(ambig_nucl, "-")

for (i in 1:length(ambigSites)){
    if (length(grep(ambigSites[i], lib_fin$seq))==0){
    }else{
        print(paste(ambigSites[i]," a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL",grep(ambigSites[i], lib_fin$seq),sep=""))
    }
}
# si aucun texte n'apparaît, c'est que tout est beau!
