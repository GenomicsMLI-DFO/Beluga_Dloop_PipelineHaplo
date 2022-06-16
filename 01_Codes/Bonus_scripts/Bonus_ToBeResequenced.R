library(readxl)
library(dplyr)

d <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "D-Loop", na = "NA"))
s <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Specimens", na = "NA"))
g <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Groupe", na = "NA"))


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
dt <- dt[dt$Nom_commun %in% "Beluga",]  # keep beluga only

# which specimens actually have an haplo because they have been sequences multiple times
rem <- which(dt$Numero_unique_specimen %in% d$Numero_unique_specimen[!is.na(d$haplotype_615)])
dt <- dt[-rem,]  # 124 specimens that are extracted once (or even multiple times) but do not have any haplotype associated with them

# nchar(d[d$Numero_unique_specimen %in% "S_20_03210", "Sequence_consensus"])

# # Prepare list of samples to resamples with info on where to find them
# ext <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Extraits_ADN_ARN", na = "NA"))
# tis <- as.data.frame(read_excel("../ACCESS/20220330_MOBELS.xlsx", sheet = "Tissus", na = "NA"))
# 
# str(ext)
# colnames(ext)[2:4] <- c("Numero_unique_specimen","Numero_unique_tissu","Numero_unique_extrait")
# ext <- subset(ext, select = c(Numero_unique_specimen,Numero_unique_tissu,Numero_unique_extrait,Numero_boite_extrait))
# 
# str(tis)
# colnames(tis)[2] <- "Numero_unique_specimen"
# tis <- subset(tis, select = c(Numero_unique_specimen,Numero_unique_tissu,No_Local_entreposage_tissu,No_congelateur_frigo,No_boite_entreposage_tissu))
# 
# samples <- left_join(dt, tis, by = c("Numero_unique_specimen"))






