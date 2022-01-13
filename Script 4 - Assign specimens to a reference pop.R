# Prepare sequences for updating haplotype libraries: long (HL, 615bp) and short (HS, 234 bp)
# Assign specimens to reference population
# 
# 
# Luca Montana
# 2021-12-17
#

getwd()
rm(list = ls())


# Libraries and function --------------------------------------------------
library(ape)
library(adegenet)
library(assignPOP)

writeGenPop <- function(gi, file.name, comment) { # from Romunov/zvau/R/writeGenPop.R
    
    if (is.list(gi)) {
        # do all genind objects have the same number of loci?
        if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
        gi.char <- gi[[1]]
        loc.names <- locNames(gi[[1]])
    } else {
        gi.char <- gi
        loc.names <- locNames(gi)
    }
    
    # Calculate the length of two alleles.
    lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
    lng <- unique(nchar(lng))
    
    stopifnot(length(lng) == 1)
    
    cat(paste(comment, "\n"), file = file.name)
    cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
    
    if (is.list(gi)) {
        pop.names <- seq_len(length(gi))
    } else {
        pop.names <- popNames(gi)
    }
    
    for (i in pop.names) {
        cat("pop\n", file = file.name, append = TRUE)
        if (is.list(gi)) {
            intm <- gi[[i]]
            loc.names <- locNames(gi[[i]])
        } else {
            intm <- gi[pop(gi) == i, drop = FALSE]
        }
        ind.names <- indNames(intm)
        intm <- genind2df(intm, sep = "")
        intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
        out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
        write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
    }
    
    return(NULL)
}


# Upload data -------------------------------------------------------------
out615 <- read.csv("Dloop_haplo615_n3283.csv")
ambiguous <- out615[which(!is.na(out615$seq)),]


connus <- read.table("Sequence_Dloop_Ref_615pb_n274_4pop.txt")
connus <- connus[-1,]
colnames(connus) <- c("id", "pop","seq")

write.csv(data2, "Dloop_haplo615_n3283.csv", row.names=F)



