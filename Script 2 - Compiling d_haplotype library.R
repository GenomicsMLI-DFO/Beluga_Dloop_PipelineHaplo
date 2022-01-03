# Prepare sequence dataset for haplotype assignment:
# Compiling the d_haplotypes library
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
data <- read.table("Sequences_Dloop_all_n3284.txt", header = T)
str(data)
colnames(data)[2] <- "seq"

# Include info on quality of sequences ------------------------------------
sequences <- toupper(data$seq)
nucleotides <- c("A","T","C","G")
ambig_nucl <- c("N","R","Y","K","M","S","W","B","D","H","V")
seq_len_exp <- 615  # valeur de la longueur attendue de la séquence

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
dna <- fasta2DNAbin("Beluga_615bp_n3284.fasta")
dna <- DNAbin2genind(dna, polyThres=0)
dna #pour information
snpPos <- locNames(dna)

# obtain info on minimal sequence
seq_start <- as.numeric(snpPos[1])  # position du premier site polymorphe
seq_stop <- as.numeric(tail(snpPos, 1))  # position du dernier site polymorphe

util <- data.frame(matrix(ncol=1, nrow=0))
colnames(util) <- c("seq_utilisable")

for (i in 1:length(sequences)){
    if (is.na(sequences[i])){
        util[i,] <- "no"
    } else if (sum(str_count(substr(sequences[i],seq_start,seq_stop),nucleotides)) == seq_stop-seq_start+1) {
        util[i,] <- "yes"  # sum no of ACTG should = minimal sequence (570 nt). Rationale: if ambiguities are present sum of A,C,T,G < 570
    }else{
        util[i,] <- "no"
    }
}
data2 <- cbind(data1, util)
write.csv(data2, "Sequences_Dloop_all_n3284.csv", row.names = F)


# Expand haplotype library ------------------------------------------------
lib <- read.csv("~/Documents/Post-Docs/IML/MOBELS/dloop/DB/dloop/librairie_139_haplotypes.csv", )
colnames(lib) <- c("hapl","seq") # if it's not already the case
table(nchar(lib$seq))  # 255 = 1; 615 = 138: there is one haplo identified on 255 nt

# New haplo library from scratches
lib <- data.frame(matrix(ncol = 2, nrow = 1))
lib[1,1] <- "HL001"
lib[1,2] <- data2$seq[data2$ID=="S_20_00600"]
colnames(lib) <- c("hapl","seq")


# detect new haplotypes
a <- length(lib$seq) # nb haplo in starting library
nlib <- data.frame(matrix(ncol=2, nrow=0))

for (i in 1:length(data2$seq)){
    if(data2$seq_utilisable[i] == "yes"){
        if (substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)){

        }else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(nlib[,2], seq_start, seq_stop)){

        }else{
            a=a+1
            if (a < 10){
                nlib <- rbind(nlib, cbind(paste("HL00",a,sep=""), data2$seq[i]))
            } else if (a < 100 & a >= 10){
                nlib <- rbind(nlib, cbind(paste("HL0",a,sep=""), data2$seq[i]))
            }else{
                nlib <- rbind(nlib,cbind(paste("HL",a,sep=""), data2$seq[i]))
            }
        }
    }else{
    }
}
colnames(nlib) <- c("hapl","seq")
lib_fin <- rbind(lib,nlib)

# Save the updated haplotype library
write.csv(lib_fin, file=paste("librairie_", length(lib_fin$hapl), "_haplotypes.csv", sep=""), row.names=F)

#data2[data2$seq==lib_fin[lib_fin$hapl=="HL104","seq"], "ID"]


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
# [1] "R a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL126"
# [1] "S a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL68"
# [1] "- a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL104"
# HL104 only found for S_20_01854
data2[data2$seq=="ACTACGTCAGTATTAAATAAACCTATTTCCAATACATTTTACTGTGACTATTGCATACCCTTATACACACACCATTAAATCTTAGTCTTTCTTTATAAATATTCATATACATATATACTATGTATTATTGTGCATTCATTTATTTTCCATACGGTCAGTTAAAGCTCGTATTAGATCTTATTAATTTTACAAATCACATAATATGCATGCTCTTACATATTATATATCAACAGTCCATTTTACCTCCATTATATACTATGGCCGCTCCATTAGATCACGAGCTTAACTACCATGCCGCGTGAAACCAGCAACCCGCTCGGCAGGGATCCCTCTTCTCGCACCGGGCCCATATCTCGTGGGGGTAGCTAATAGTGGTCTTCACAAGACATCTGGTTCTTACTTCAGGACCATTCCAGCTTAAAATCGCCCACTCGTTCCCCTTAAATAAGACATCTCGATGGACTAATGACTAATCAGCCCATGCTCACACATAACTGAGATTTCATACATTTGGTATTTTTTATTTTTGGGGGGGGGCCTGCACCGACTCAGCTATGGCCTTAGAAAGGCCCTGTCACAGTCAGATAAATTGTAGCTGGACCTGTGTGTATTTTT",]
# HL68 also found in S_20_01102-2, S_20_01436, S_20_02110, S_20_02167, S_20_02568, S_20_03122, S_20_03198, S_20_03652, S_20_03679 and with C instead of S
data2[data2$seq=="ACTACGTCAGTATTAAATAAACCTATTTCCAATACATTTTACTGTGACTATTGCATACCCTTATACACACACCATTAAATCTTAGTCTTTCTTTATAAATATTCATATACATATATACTATGTATTATTGTGCATTCATTTATTTTCCATACGGTCAGTTAAAGCTCGTATTAGATATTATTAATTTTACAAATCACATAATATGCATGCTCTTACATATTATATATCAACAGTCCATTTTACCTCCATTATATACTATGGCCGCTCCATTAGATCACGAGCTTAACTACCATGCCGCGTGAAACCAGCAACCCGCTCGGCAGGGATCCCTCTTCTCGCACCGGGCCCATATCTCGTGGGGGTAGCTAATAGTGGTCTTCACAAGACATCTGGTTCTTACTTCAGGACCATTCCAACTTAAAATCGCCCACTCGTTCCCCTTAAATAAGACATCTCGATGGACTAATGACTAATCAGCCCATGCTCACACATAACTGAGATTTCATACATTTGGTATTTTTTATTTTTGGGGGGGGGCCTGCACCGACTCAGCTATGGCCTTAGAAAGGCCCTGTCACAGTCAGATAAATTGTAGCTGGACCTGTGTGTATTTTT",]
# HL126 also found for S_20_03546 with A instead of R

# Change sequences in lib haplo with ambiguitis (outside minimal sequence) with 'clean' sequence
lib_fin[lib_fin$hapl=="HL068", "seq"] <- data2$seq[data2$ID=="S_20_01102-2"]
lib_fin[lib_fin$hapl=="HL126", "seq"] <- data2$seq[data2$ID=="S_20_03546"]
write.csv(lib_fin, file=paste("librairie_", length(lib_fin$hapl), "_haplotypes.csv", sep=""), row.names=F)
