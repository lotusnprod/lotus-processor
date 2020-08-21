# title: "sampleR"

# loading functions
source("functions.R")
source("paths.R")

# loading files
print(x = "loading db, if running fullmode, this may take a while")

## fullDB
openDb <- read_delim(
  file = gzfile(pathDataInterimTablesAnalysedOpenDbTriplets),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

referenceTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB_doi <- openDb %>%
  filter(referenceType == "doi") %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB_original <- openDb %>%
  filter(referenceType == "original")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

# set.seed(seed = 42,
#          kind = "Mersenne-Twister",
#          normal.kind = "Inversion") # 0 entries
# sampleONPDB_publishingDetails <- openDb %>%
#   filter(referenceType == "publishingDetails") %>%
#   sample_n(30) %>%
#   mutate(
#     curator = NA,
#     validated = NA,
#     comments = NA
#   )

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB_pubmed <- openDb %>%
  filter(referenceType == "pubmed") %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB_split <- openDb %>%
  filter(referenceType == "split")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB_title <- openDb %>%
  filter(referenceType == "title")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

sampleONPDB <- bind_rows(
  sampleONPDB_doi,
  sampleONPDB_original,
  # sampleONPDB_publishingDetails,
  sampleONPDB_pubmed,
  sampleONPDB_split,
  sampleONPDB_title
) %>%
  mutate_all(as.character)


set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB <- sampleONPDB[sample(nrow(sampleONPDB)), ]

sampleONPDB[1:50, "curator"] <- "AR"

sampleONPDB[51:100, "curator"] <- "JB"

sampleONPDB[101:150, "curator"] <- "PMA"

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleKnapsack <- openDb %>%
  filter(database == "kna_1") %>%
  sample_n(150) %>%
  mutate(
    curator = sample(c("AR", "JB", "PMA"),
                     size = nrow(.),
                     replace = TRUE),
    validated = NA,
    comments = NA
  )

openDbFull <- left_join(openDb, referenceTableFull)

goldenSet <- openDbFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmid) |
      !is.na(referenceCleanedPmcid)
  ) %>%
  filter(
    referenceCleaned_score_crossref == 1 |
      referenceCleaned_score_distance <= 10 |
      # here is a discussion about | or &
      referenceCleaned_score_titleOrganism == 1
  ) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
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
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

platinumSet <- openDbFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmid) |
      !is.na(referenceCleanedPmcid)
  ) %>%
  filter(referenceCleaned_score_crossref == 1 |
           referenceCleaned_score_distance <= 10) %>%
  filter(referenceCleaned_score_titleOrganism == 1) %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
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
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid
  )

set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleWD <- goldenSet %>%
  sample_n(500) %>%
  select(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
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
  file = gzfile(
    description = pathDataInterimTablesAnalysedGold,
    compression = 9,
    encoding = "UTF-8"
  ),
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

## platinumSet
write.table(
  x = platinumSet,
  file = pathDataInterimTablesAnalysedPlatinum,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
