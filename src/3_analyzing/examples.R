source("r/log_debug.R")
log_debug("This script helps finding nice examples")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("loading ...")
log_debug("... validated db, if running fullmode, this may take a while")
openDb <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedPlatinum),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

structureSearch_1 <- openDb %>%
  filter(structureType == "nominal") %>%
  distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch_2 <- openDb %>%
  filter(structureType == "smiles") %>%
  distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch_3 <- openDb %>%
  filter(structureType == "inchi") %>%
  distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch <-
  rbind(structureSearch_1, structureSearch_2, structureSearch_3) %>%
  group_by(structureCleanedInchikey) %>%
  add_count() %>%
  arrange(desc(n))

saltSearch <- structureSearch_3 %>%
  distinct(structureValue,
    .keep_all = TRUE
  ) %>%
  group_by(structureCleanedInchikey) %>%
  add_count() %>%
  filter(grepl(pattern = "\\.", x = structureValue)) %>%
  arrange(desc(n))

maybeHit_salt <- openDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismCleaned)) %>%
  filter(structureCleanedInchikey == "KRKNYBCHXYNGOX-UHFFFAOYSA-N") %>%
  distinct(structureValue)

maybeHit_str <- openDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismCleaned)) %>%
  filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N")

hitNames_str <- maybeHit_str %>%
  filter(structureType == "nominal") %>%
  distinct(structureValue)

hitSmiles_str <- maybeHit_str %>%
  filter(structureType == "smiles") %>%
  distinct(structureValue)

hitInchi_str <- maybeHit_str %>%
  filter(structureType == "inchi") %>%
  distinct(structureValue)

organismSearch <- openDb %>%
  distinct(organismType,
    organismValue,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  group_by(organismCleaned) %>%
  add_count() %>%
  arrange(desc(n)) %>%
  filter(
    !grepl(pattern = "Streptomyces", x = organismCleaned) &
      !grepl(pattern = "Aspergillus", x = organismCleaned) &
      !grepl(pattern = "Fusarium", x = organismCleaned) &
      !grepl(pattern = ".*ae", x = organismCleaned)
  )

maybeHit_org <- openDb %>%
  filter(organismCleaned == "Oryza sativa")

hitNames_org <- maybeHit_org %>%
  distinct(
    organismType,
    organismValue
  )

reference <- openDb %>%
  filter(!is.na(referenceCleanedDoi))

referenceSearch <- reference %>%
  distinct(referenceValue,
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

maybeHit_ref <- openDb %>%
  filter(referenceCleanedDoi == "10.1021/np0600595")

hitNames_org <- maybeHit_org %>%
  distinct(
    organismType,
    organismValue
  )

doubleTest <- openDb %>%
  filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") %>%
  distinct(
    organismType,
    organismValue, organismCleaned
  ) %>%
  group_by(organismCleaned) %>%
  count() %>%
  arrange(desc(n))

pairTest <- openDb %>%
  filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") %>%
  filter(organismCleaned == "Crataegus monogyna") %>%
  distinct(
    organismType,
    organismValue, structureValue, organismCleaned
  )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
