library(dplyr)
library(readr)

source("paths.R")

inhouseDb <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

structure <- inhouseDb %>%
  filter(!is.na(inchikeySanitized))

structureSearch_1 <- structure %>%
  filter(is.na(structureOriginalInchi)) %>%
  filter(is.na(structureOriginalSmiles)) %>%
  filter(!is.na(structureOriginalNominal)) %>%
  distinct(structureOriginalNominal,
           inchikeySanitized, .keep_all = TRUE)

structureSearch_2 <- structure %>%
  filter(is.na(structureOriginalInchi)) %>%
  filter(!is.na(structureOriginalSmiles)) %>%
  distinct(structureOriginalSmiles,
           inchikeySanitized, .keep_all = TRUE)

structureSearch_3 <- structure %>%
  filter(!is.na(structureOriginalInchi)) %>%
  distinct(structureOriginalInchi,
           inchikeySanitized, .keep_all = TRUE)

structureSearch <-
  rbind(structureSearch_1, structureSearch_2, structureSearch_3)

structureSearch <- structureSearch %>%
  group_by(inchikeySanitized) %>%
  add_count() %>%
  arrange(desc(n))

saltSearch <- structureSearch_3 %>%
  distinct(structureOriginalInchi, .keep_all = TRUE) %>%
  group_by(inchikeySanitized) %>%
  add_count() %>%
  filter(grepl(pattern = "\\.", x = structureOriginalInchi)) %>%
  arrange(desc(n))

maybeHit_salt <- inhouseDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismLowestTaxon)) %>%
  filter(inchikeySanitized == "KRKNYBCHXYNGOX-UHFFFAOYSA-N") %>%
  distinct(structureOriginalInchi)

maybeHit_str <- inhouseDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismLowestTaxon)) %>%
  filter(inchikeySanitized == "XMGQYMWWDOXHJM-UHFFFAOYSA-N")

hitNames_str <- maybeHit_str %>%
  distinct(structureOriginalNominal)

hitSmiles_str <- maybeHit_str %>%
  distinct(structureOriginalSmiles)

hitInchi_str <- maybeHit_str %>%
  distinct(structureOriginalInchi)



organism <- inhouseDb %>%
  filter(!is.na(organismLowestTaxon))

organismSearch <- organism %>%
  distinct(organismOriginal,
           organismLowestTaxon,
           .keep_all = TRUE) %>%
  group_by(organismLowestTaxon) %>%
  add_count() %>%
  arrange(desc(n)) %>%
  filter(
    !grepl(pattern = "Streptomyces", x = organismLowestTaxon) &
      !grepl(pattern = "Aspergillus", x = organismLowestTaxon) &
      !grepl(pattern = "Fusarium", x = organismLowestTaxon) &
      !grepl(pattern = ".*ae", x = organismLowestTaxon)
  )

maybeHit_org <- inhouseDb %>%
  filter(organismLowestTaxon == "Oryza sativa")

hitNames_org <- maybeHit_org %>%
  distinct(organismOriginal)



reference <- inhouseDb %>%
  filter(!is.na(referenceCleanedDoi))

referenceSearch <- reference %>%
  distinct(
    referenceOriginalDoi,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    referenceCleanedDoi,
    .keep_all = TRUE
  ) %>%
  group_by(referenceCleanedDoi) %>%
  add_count() %>%
  arrange(desc(n)) %>%
  filter(
    !grepl(pattern = "Khimiya", x = referenceCleanedTitle) &
      !grepl(pattern = "Flavone and flavonol glycosides", x = referenceCleanedTitle) &
      !grepl(pattern = "The Handbook of Natural Flavonoids", x = referenceCleanedTitle) &
      n <= 10 &
      n >= 3
  )

maybeHit_ref <- inhouseDb %>%
  filter(referenceCleanedDoi == "10.1021/np0600595")

hitNames_org <- maybeHit_org %>%
  distinct(organismOriginal)
