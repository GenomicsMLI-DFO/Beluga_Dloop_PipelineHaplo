# Prepare sequence dataset for haplotype assignment:
# Compiling the d_haplotypes library
# 
# Benjamin Hornoy
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


# Include info on quality of sequences ------------------------------------
sequences <- toupper(data$Sequence)
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
seq_start <- as.numeric(snpPos[1])  # position du premier site polymorphe
seq_stop <- as.numeric(tail(snpPos, 1))  # position du dernier site polymorphe

util <- data.frame(matrix(ncol=1, nrow=0))
colnames(util) <- c("seq_utilisable")

for (i in 1:length(sequences)){
    if (is.na(sequences[i])){
        util[i,] <- "0"
    } else if (sum(str_count(substr(sequences[i],seq_start,seq_stop),nucleotides)) == seq_stop-seq_start+1) {
        util[i,] <- "1"			
    }else{
        util[i,] <- "0"
    }
}
data2 <- cbind(data1, util)



#### CHECK THIS OUT ####
# Expand haplotype library ------------------------------------------------
# lib <- read.csv("librairie_107_haplotypes.csv")
# colnames(lib) <- c("hapl","seq") # si ce n'est pas déjà le cas
# 
# # charger les infos sur la séquence minimale
# seq_start <- 22 #position du premier site polymorphe
# seq_stop <- 605 #position du dernier site polymorphe
# 
# # détecter les nouveaux haplotypes
# a=length(lib$seq) # nb haplo dans la librairie initiale
# library <- data.frame(matrix(ncol=2, nrow=0))
# 
# for (i in 1:length(data2$seq)){
#     if(data2$seq_utilisable[i] == "oui"){
#         if (substr(data2$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)){
#             
#         }else if(substr(data2$seq[i], seq_start, seq_stop) %in% substr(library[,2], seq_start, seq_stop)){
#             
#         }else{
#             a=a+1
#             if (a < 100){
#                 library <- rbind(library, cbind(paste("HL0",a,sep=""), data2$seq[i]))
#             }else{
#                 library <- rbind(library,cbind(paste("HL",a,sep=""), data2$seq[i]))
#             }
#         }
#     }else{		
#     }
# }
# 
# colnames(library) <- c("hapl","seq")
# lib_fin <- rbind(lib,library)
# 
# # sauver la librairie enrichie
# write.csv(lib_fin, file=paste("librairie_", length(lib_fin$hapl), "_haplotypes.csv", sep=""), row.names=F)
# 
# ### ATTENTION! La "seq" peut contenir des sites ambigus ou manquants (hors de la séquence minimale), ce qui pourra poser pb dans le futur si ces sites sont détectés comme polymorphes dans le futur. Il faut donc automatiser ou s'assurer à la mitaine que la séquence mise dans la librairie ne contient pas d'ambigus, en espérant qu'aucun haplotype ne soit représenté que par un seul individu, qui a des ambigus...
# 
# # vérifier si des ambigus sont présents dans la nouvelle librairie (si oui, voir ci-dessus...)
# ambigSites <- c(ambig_nucl, "-")
# 
# for (i in 1:length(ambigSites)){
#     if (length(grep(ambigSites[i], lib_fin$seq))==0){
#     }else{
#         print(paste(ambigSites[i]," a été trouvé dans la/les séquence(s) représentative(s) du/des haplotype(s): HL",grep(ambigSites[i], lib_fin$seq),sep=""))
#     }
# }
# # si aucun texte n'apparaît, c'est que tout est beau!
