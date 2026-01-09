library(dplyr)
library(stringr)
library(tibble)

verif_consensus_sequences <- function(dt,
                                      seq_col = "Sequence_consensus",
                                      min_length = 200) {
  nt <- c("A", "T", "C", "G")
  ambiguous <- c("N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V")

  # Séquences en majuscule (sans modifier l'objet original)
  dt <- dt %>% mutate(!!seq_col := toupper(.data[[seq_col]]))

  # Longueur
  dt <- dt %>% mutate(N_nucl = nchar(.data[[seq_col]]))

  # Audit des positions ambiguës
  find_ambiguous_pos <- function(seq) {
    if (is.na(seq)) return(NA)
    v <- strsplit(seq, "")[[1]]
    amb_pos <- which(v %in% ambiguous)
    if (length(amb_pos) == 0) return(NA)
    paste(amb_pos, collapse = ",")
  }
  dt <- dt %>% mutate(Ambiguous_positions = sapply(.data[[seq_col]], find_ambiguous_pos))

  # Vérifications
  pb_dash <- dt %>% filter(str_detect(.data[[seq_col]], "-")) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = "Contient '-'")

  pb_colon <- dt %>% filter(str_detect(.data[[seq_col]], ":")) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = "Contient ':'")

  pb_short <- dt %>%
    filter(N_nucl < min_length) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = paste0("< ", min_length, " nucléotides"))

  pb_na <- dt %>% filter(is.na(.data[[seq_col]])) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = "NA Sequence")

  base_allowed <- c(nt, ambiguous)
  pb_nonstandard <- dt %>%
    filter(str_detect(.data[[seq_col]], paste0("[^", paste(base_allowed, collapse=""), "]"))) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = "Nucléotides non standards (hors ATCGNRYKMSWBDHV)")

  pb_ambiguity <- dt %>%
    filter(str_detect(.data[[seq_col]], paste(ambiguous[ambiguous != "N"], collapse="|"))) %>%
    select(Numero_unique_specimen, Numero_unique_extrait, ID, Ambiguous_positions) %>%
    mutate(Anomalie = "Base(s) ambiguë(s) (RYKMSWBDHV)")

  # Résumé
  res <- bind_rows(pb_dash, pb_colon, pb_short, pb_na, pb_nonstandard, pb_ambiguity)

  if (nrow(res) == 0) {
    message("Tout est OK dans la colonne '", seq_col, "' : aucune modification requise.")
    return(invisible(NULL))
  } else {
    message("Des modifications sont requises sur certaines lignes :")
    return(res)
  }
}
