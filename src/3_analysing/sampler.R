# title: "sampleR"

# loading functions
source("functions.R")
source("paths.R")

# loading files
print(x = "loading db, if running fullmode, this may take a while")

## fullDB
openDbTriplets <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedOpenDbTriplets),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

set.seed(42)
sampleONPDB <- openDbTriplets %>%
  sample_n(150) %>%
  mutate(
    curator = sample(c("AR", "JB", "PMA"),
                     size = nrow(.),
                     replace = TRUE),
    validated = NA,
    comments = NA
  )

set.seed(42)
sampleKnapsack <- openDbTriplets %>%
  filter(database == "kna_1") %>%
  sample_n(150) %>%
  mutate(
    curator = sample(c("AR", "JB", "PMA"),
                     size = nrow(.),
                     replace = TRUE),
    validated = NA,
    comments = NA
  )

goldenSet <- openDbTriplets %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmid) |
      !is.na(referenceCleanedPmcid)
  ) %>%
  filter(referenceCleaned_score_crossref == 100 |
           referenceCleaned_score_distance <= 10) %>%
  filter(referenceCleaned_score_titleOrganism == 1) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonId,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  ) %>%
  select(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonId,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

set.seed(42)
sampleWD <- goldenSet %>%
  sample_n(500) %>%
  select(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonId,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

#exporting
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesAnalysed),
  dir.create(pathDataInterimTablesAnalysed),
  FALSE
)

## sampleONPDB
write.table(
  x = sampleONPDB,
  file = pathDataInterimTablesAnalysedSampleAllONPDB,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## sampleKnapsack
write.table(
  x = sampleKnapsack,
  file = pathDataInterimTablesAnalysedSampleKnapsack,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## goldenSet
write.table(
  x = goldenSet,
  file = pathDataInterimTablesAnalysedGold,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## sampleWD
write.table(
  x = sampleWD,
  file = pathDataInterimTablesAnalysedSampleGoldONPDB,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
