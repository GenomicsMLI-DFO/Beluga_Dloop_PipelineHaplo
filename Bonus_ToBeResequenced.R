library(readxl)
library(dplyr)

d <- read_excel("~/MOBELS/DB/ACCESS/20220329_MOBELS.xlsx", sheet = "D-Loop", na = "NA")
s <- read_excel("~/MOBELS/DB/ACCESS/20220329_MOBELS.xlsx", sheet = "Specimens", na = "NA")
g <- read_excel("~/MOBELS/DB/ACCESS/20220329_MOBELS.xlsx", sheet = "Groupe", na = "NA")


str(s)
s <- s[, names(s) %in% c("Numero_unique_specimen","Nom_commun","")]

colnames(d)[2] <- "Numero_unique_extrait"
colnames(d)
d <- subset(d, select = c(Numero_unique_specimen,Numero_unique_extrait,Qualite_sequence,Sequence_utilisable_234,Sequence_utilisable_615,Sequence_consensus,N_nucl_615,N_ambig_615,
                          N_manquants_615,haplotype_615))
colnames(s)
s <- subset(s, select = c(Numero_unique_specimen,Nom_commun))
colnames(g)
g <- subset(g, select = c(Numero_unique_reception_groupe,Region_echantillonnage,Lieu_echantillonnage,Annee_echantillonnage:Jour_echantillonnage))

dt <- left_join(subset(d, select = -Sequence_consensus), s, by = "Numero_unique_specimen")
dt <- left_join(dt, g, by = c("Numero_unique_specimen"="Numero_unique_reception_groupe"))
dt <- dt[is.na(dt$haplotype_615),]
table(dt$Region_echantillonnage, useNA = "ifany")
dt <- dt[!is.na(dt$Region_echantillonnage),]
dt <- dt[!is.na(dt$Mois_echantillonnage),]

# which specimens actually have an haplo because they have been sequences multiple times
rem <- which(dt$Numero_unique_specimen %in% d$Numero_unique_specimen[!is.na(d$haplotype_615)])
dt <- dt[-rem,]
