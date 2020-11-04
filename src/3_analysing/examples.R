cat("This script helps finding nice examples \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions.R")

cat("loading ... \n")
cat("... validated db, if running fullmode, this may take a while \n")
openDb <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedPlatinum),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

structureSearch_1 <- openDb %>%
  filter(structureType == "nominal") %>%
  distinct(structureValue,
    structureCleanedInchikey3D,
    .keep_all = TRUE
  )

structureSearch_2 <- openDb %>%
  filter(structureType == "smiles") %>%
  distinct(structureValue,
    structureCleanedInchikey3D,
    .keep_all = TRUE
  )

structureSearch_3 <- openDb %>%
  filter(structureType == "inchi") %>%
  distinct(structureValue,
    structureCleanedInchikey3D,
    .keep_all = TRUE
  )

structureSearch <-
  rbind(structureSearch_1, structureSearch_2, structureSearch_3) %>%
  group_by(structureCleanedInchikey3D) %>%
  add_count() %>%
  arrange(desc(n))

saltSearch <- structureSearch_3 %>%
  distinct(structureValue,
    .keep_all = TRUE
  ) %>%
  group_by(structureCleanedInchikey3D) %>%
  add_count() %>%
  filter(grepl(pattern = "\\.", x = structureValue)) %>%
  arrange(desc(n))

maybeHit_salt <- openDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismCleaned)) %>%
  filter(structureCleanedInchikey3D == "KRKNYBCHXYNGOX-UHFFFAOYSA-N") %>%
  distinct(structureValue)

maybeHit_str <- openDb %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(!is.na(organismCleaned)) %>%
  filter(structureCleanedInchikey3D == "OVSQVDMCBVZWGM-DTGCRPNFSA-N")

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
  distinct(organismOriginal,
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
  distinct(organismOriginal)

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
  distinct(organismOriginal)

doubleTest <- openDb %>%
  filter(structureCleanedInchikey3D == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") %>%
  distinct(organismOriginal, organismCleaned) %>%
  group_by(organismCleaned) %>%
  count() %>%
  arrange(desc(n))

pairTest <- openDb %>%
  filter(structureCleanedInchikey3D == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") %>%
  filter(organismCleaned == "Crataegus monogyna") %>%
  distinct(organismOriginal, structureValue, organismCleaned)