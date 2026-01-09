# Info --------------------------------------------------------------------
#
# Authors: Benjamin Hornoy, Luca Montana
# Updated by: Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch
# Location: Institut Maurice Lamontagne
# Date: 2021-2026
#
# Overview: Assign haplotypes to specimens

# 0. Housekeeping ---------------------------------------------------------

# Verify if you're in the right directory
#getwd()

# Libraries
#library(readxl)
library(dplyr)
library(tidyverse)

# Functions
"%nin%" <- Negate("%in%")


# 1. Data -----------------------------------------------------------------

projet <- "Monitorage_20260108"

if (!dir.exists(file.path("02_Results/02_ACCESS/", projet)))
  dir.create(file.path("02_Results/02_ACCESS/", projet))

metadata <- read_csv(paste0("00_Data/02_dloop_clean/", projet, "/", projet, "_metadata.csv"))

# data includes info on quality of sequences (columns N.nucl, N.ATCG, N.ambog, N.manquants) as well as if sequences is usable
# all made in 2a_HaploLibrary_234.R and 2b_HaploLibrary_615.R

s570.red <- Biostrings::readDNAStringSet(paste0("00_Data/02_dloop_clean/", projet, "/", projet, "_570bp.fasta"))
s570.red

## 2.1. Upload haplotype libraries ------------------------------------------

list.files("02_Results/00_libraries/")

lib570.fasta <- Biostrings::readDNAStringSet('02_Results/00_libraries/librairies_148_haplotype570.fasta')

lib570 <- data.frame(hapl_570 = names(lib570.fasta),
                     seq = as.character(lib570.fasta))


# 2. Assign haplotype to each individual ----------------------------------
s570.df <- data.frame(ID = names(s570.red),
                      seq = as.character(s570.red) |> str_remove_all("-")) |>
  dplyr::left_join(lib570) |>
  #dplyr::select(-seq)|>
  dplyr::mutate(hapl_570 = ifelse(is.na(hapl_570), "NewHap", hapl_570))


# Is there new haplotypes?
s570.df |> dplyr::filter(hapl_570 == "NewHap")


complete.metadata.df <-  metadata |> mutate(ID = as.character(ID))|> left_join(s570.df |> dplyr::select(-seq), by = c("ID"))

library(ggplot2)
complete.metadata.df |> dplyr::group_by(hapl_570, Region_echantillonnage, Lieu_echantillonnage) |>
  summarise(N = n()) |>
  ggplot(aes(x= Lieu_echantillonnage, y = hapl_570, fill = N)) +
  geom_bin2d() +
  facet_grid(.~Region_echantillonnage, space = "free", scale = "free") +
  scale_fill_viridis_c(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


complete.metadata.df |> dplyr::filter(!is.na(hapl_570)) |>
  group_by(Numero_unique_specimen, Numero_unique_extrait) |>
  summarise(N = length(hapl_570),
            Ndiff = length(unique(hapl_570)),
            hapl = paste(unique(hapl_570), collapse = ",")) |>
  arrange(desc(Ndiff)) |> View()

readr::write_csv(complete.metadata.df, file.path("02_Results/02_ACCESS/", projet, paste0(projet, "_haplo.csv")))
