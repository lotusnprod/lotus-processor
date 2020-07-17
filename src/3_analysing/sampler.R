# title: "sampleR"

# loading functions
source("functions.R")
source("paths.R")

# loading files
print(x = "loading db, if running fullmode, this may take a while")

## fullDB
inhouseDb <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(referenceCleanedTranslationScore = as.integer(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceOriginalExternal)) %>%
  arrange(desc(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceCleanedDoi)) %>% #very important to keep references
  data.frame()

openDbTriplets <- inhouseDb %>%
  filter(!is.na(inchikeySanitized)) %>%
  filter(!is.na(organismLowestTaxon)) %>%
  filter(
    !is.na(referenceCleanedTitle) |
      !is.na(referenceCleanedDoi) |
      referenceOriginalExternal == "PLANTCYC" |
      referenceOriginalExternal == "Korean_Traditional_Knowledge_Portal_20140617"
  ) %>%
  select(
    database,
    organismOriginal,
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    organismLowestTaxon,
    inchikeySanitized,
    smilesSanitized,
    referenceCleanedTitle,
    referenceCleanedDoi,
    referenceCleanedTranslationScore
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

set.seed(42)
sampleWD <- openDbTriplets %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(referenceCleanedTranslationScore == 100) %>%
  sample_n(1000) %>%
  select(
    organismLowestTaxon,
    inchikeySanitized,
    smilesSanitized,
    referenceCleanedTitle,
    referenceCleanedDoi
  )

#exporting
print(x = "exporting, may take a while if running full mode")
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


## sampleWD
write.table(
  x = sampleWD,
  file = pathDataInterimTablesAnalysedSampleGoldONPDB,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
