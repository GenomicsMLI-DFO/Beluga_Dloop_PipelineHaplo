# Prepare sequence dataset for haplotype assignment:
# Compiling the d_haplotypes library
# 
# Benjamin Hornoy
#
#

getwd()
rm(list = ls())

# Libraries and function --------------------------------------------------


# Data --------------------------------------------------------------------
data <- read.csv("Sequences_Dloop_all_n3284.csv")
str(data)
# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in Script 2 - Compiling d_haplotype library.R


# Assign haplotype to each individuals ------------------------------------
lib <- read.csv('librairie_136_haplotypes.csv') # most recent haplotype library
colnames(lib) <- c("hapl","seq")  # if it's not already the case

# Upload info on minimal sequence
min_seq <- read.csv("polymorphismes_et_seq_minimale_n3102.csv", stringsAsFactors = F)  # table made in script 1
seq_start <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale, " - "))[1])  # start of minimal sequence
seq_stop <- as.numeric(unlist(strsplit(min_seq$Bornes.sequence.minimale, " - "))[2])  # end of minimal sequence

# Assign haplotype to each individual
hapind <- data.frame(matrix(ncol=1, nrow=0))
colnames(hapind) <- c("haplotype")

sequences <- toupper(data$seq)
for (i in 1:length(sequences)){
    if(data$seq_utilisable[i] == "yes"){
        if (substr(data$seq[i], seq_start, seq_stop) %in% substr(lib$seq, seq_start, seq_stop)){
            hapind[i,1] <- lib$hapl[which(substr(lib$seq, seq_start, seq_stop) == substr(data$seq[i], seq_start, seq_stop))]
        }else{
            # shows if haplotype is new relative to library
            hapind[i,1] <- "unknown haplotype"
        }
    }else{
        # generates NAs for unusable sequences, to avoid confusion
        hapind[i,1] <- "NA"
    }
}

table(hapind)
data2 <- cbind(data, hapind)

write.csv(data2, "Dloop_haplo_n3284.csv", row.names=F)

