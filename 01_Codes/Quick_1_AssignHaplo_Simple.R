# Info --------------------------------------------------------------------
#
# Authors: Benjamin Hornoy, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch
# Location: Institut Maurice Lamontagne
# Date: 2021-12-17
#
# Overview: Assign haplotypes to specimens
#
#


# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
#getwd()

# Libraries
#library(readxl)
library(dplyr)

# Functions
"%nin%" <- Negate("%in%")


# 1. Data -----------------------------------------------------------------

projet <- "BelugaMonitorage_Part1_20240503"

metadata <- read.table(paste0("00_Data/02_dloop_clean/", projet, "_metadata.txt"), header = T)
str(metadata)


# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in 2a_HaploLibrary_234.R and 2b_HaploLibrary_615.R

list.files("00_Data/01_fasta/", pattern = "onlyATGC") |> str_subset(projet)

s234.red <- Biostrings::readDNAStringSet("00_Data/01_fasta/BelugaMonitorage_Part1_20240503_234bp_onlyATGC_n267.fasta")

s570.red <- Biostrings::readDNAStringSet("00_Data/01_fasta/BelugaMonitorage_Part1_20240503_570bp_onlyATGC_n266.fasta")

s234.red
s570.red

## 2.1. Upload haplotype libraries ------------------------------------------

lib234 <- read.csv('02_Results/00_libraries/librairie_54_haplotypes234.csv')  # most recent haplotype library
colnames(lib234) <- c("hapl_234","seq")  # if it's not already the case
#lib615 <- read.csv('02_Results/00_libraries/librairie_143_haplotypes615.csv')  # most recent haplotype library
#colnames(lib615) <- c("hapl_615","seq")  # if it's not already the case
lib570 <- read.csv('02_Results/00_libraries/librairie_143_haplotypes570.csv')  # most recent haplotype library
colnames(lib570) <- c("hapl_570","seq")


lib234.fasta <- Biostrings::DNAStringSet(lib234$seq)
names(lib234.fasta) <- lib234$hapl_234

lib570.fasta <- Biostrings::DNAStringSet(lib570$seq)
names(lib570.fasta) <- lib570$hapl_570



# 2. Assign haplotype to each individual ----------------------------------

s234.df <- data.frame(Numero_unique_specimen = names(s234.red),
                      seq = as.character(s234.red) |> str_remove_all("-")) |>
                      dplyr::left_join(lib234) |>
                      #dplyr::select(-seq) |>
                      dplyr::mutate(hapl_234 = ifelse(is.na(hapl_234), "NewHap", hapl_234))

#s615.df <- data.frame(Numero_unique_specimen = names(s615.red),
#                      seq = as.character(s615.red)) |> dplyr::left_join(lib615) |>
#  dplyr::select(-seq)|>
#  dplyr::mutate(hapl_615 = ifelse(is.na(hapl_615), "NewHap", hapl_615))

s234.new <- s234.df |> dplyr::filter(hapl_234 == "NewHap") |> pull(seq)
duplicated(s234.new) |> table()

lib234.fasta[c("HS055", "HS056")] <- s234.new
lib234.fasta

lib234.df <- data.frame(hapl_234 = names(lib234.fasta),
                        seq = as.character(lib234.fasta))


Biostrings::writeXStringSet(lib234.fasta, "02_Results/00_libraries/librairies_56_haplotype234.fasta")

# Ajout

s234.df <- data.frame(Numero_unique_specimen = names(s234.red),
                      seq = as.character(s234.red) |> str_remove_all("-")) |>
  dplyr::left_join(lib234.df) |>
  #dplyr::select(-seq) |>
  dplyr::mutate(hapl_234 = ifelse(is.na(hapl_234), "NewHap", hapl_234))


s234.df |> dplyr::filter(hapl_234 == "NewHap")


s570.df <- data.frame(Numero_unique_specimen = names(s570.red),
                      seq = as.character(s570.red) |> str_remove_all("-")) |>
  dplyr::left_join(lib570) |>
  #dplyr::select(-seq)|>
  dplyr::mutate(hapl_570 = ifelse(is.na(hapl_570), "NewHap", hapl_570))




s570.new <- s570.df |> dplyr::filter(hapl_570 == "NewHap") |> pull(seq)
duplicated(s570.new) |> table()

lib570.fasta[c("HL144", "HL145", "HL146")] <- s570.new
lib570.fasta

lib570.df <- data.frame(hapl_570 = names(lib570.fasta),
                        seq = as.character(lib570.fasta))


Biostrings::writeXStringSet(lib570.fasta, "02_Results/00_libraries/librairies_146_haplotype570.fasta")

s570.df <- data.frame(Numero_unique_specimen = names(s570.red),
                      seq = as.character(s570.red) |> str_remove_all("-")) |>
  dplyr::left_join(lib570.df) |>
  #dplyr::select(-seq) |>
  dplyr::mutate(hapl_570 = ifelse(is.na(hapl_570), "NewHap", hapl_570))


s570.df |> dplyr::filter(hapl_570 == "NewHap")



complete.df <- s234.df |>  dplyr::select(-seq) |> dplyr::full_join(s570.df |> dplyr::select(-seq))


complete.df |> head()
complete.df |> View()
metadata |> head()
complete.metadata.df <-  metadata |> left_join(complete.df, by = c("ID" = "Numero_unique_specimen")) |> dplyr::select(-Sequence)

library(ggplot2)
complete.metadata.df |> dplyr::group_by(hapl_570, Region_echantillonnage) |>
  summarise(N = n()) |>
  ggplot(aes(x= Region_echantillonnage, y = hapl_570, fill = N)) +
  geom_bin2d() +
  scale_fill_viridis_c(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


complete.metadata.df |> dplyr::filter(hapl_570 == "NewHap")

readr::write_csv(complete.metadata.df, "Beluga_haplo_20240506.csv")

table(complete.df$hapl_234)
table(complete.df$hapl_615)

nrow(s615.df)
nrow(s570.df)
nrow(s234.df)
