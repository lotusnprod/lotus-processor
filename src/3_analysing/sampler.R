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
) %>%
  mutate(referenceCleanedTranslationScore = as.integer(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceOriginal_external)) %>%
  arrange(desc(referenceCleanedTranslationScore)) %>%
  arrange(desc(referenceCleanedDoi)) %>% #very important to keep references
  data.frame()

openDbTripletsClean <- openDbTriplets %>%
  filter(!is.na(inchikeySanitized)) %>%
  filter(!is.na(organismLowestTaxon)) %>%
  filter(
    !is.na(referenceCleanedTitle) |
      !is.na(referenceCleanedDoi) |
      referenceOriginal_external == "PLANTCYC" |
      referenceOriginal_external == "Korean_Traditional_Knowledge_Portal_20140617"
  ) %>%
  select(
    database,
    organismOriginal,
    structureOriginal_inchi,
    structureOriginal_smiles,
    structureOriginal_nominal,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_publishingDetails,
    referenceOriginal_title,
    referenceOriginal_split,
    organismLowestTaxon,
    organismDbTaxo,
    organismTaxonId,
    inchikeySanitized,
    inchiSanitized,
    smilesSanitized,
    referenceCleanedTitle,
    referenceCleanedDoi,
    referenceCleanedTranslationScore
  ) %>%
  distinct(
    inchikeySanitized,
    organismLowestTaxon,
    referenceOriginal_external,
    referenceCleanedDoi,
    referenceCleanedTitle,
    .keep_all = TRUE
  )

set.seed(42)
sampleONPDB <- openDbTripletsClean %>%
  sample_n(150) %>%
  mutate(
    curator = sample(c("AR", "JB", "PMA"),
                     size = nrow(.),
                     replace = TRUE),
    validated = NA,
    comments = NA
  )

set.seed(42)
sampleKnapsack <- openDbTripletsClean %>%
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
sampleWD <- openDbTripletsClean %>%
  filter(!is.na(referenceCleanedDoi)) %>%
  filter(referenceCleanedTranslationScore == 100) %>%
  sample_n(1000) %>%
  select(
    organismLowestTaxon,
    organismDbTaxo,
    organismTaxonId,
    inchikeySanitized,
    inchiSanitized,
    smilesSanitized,
    referenceCleanedTitle,
    referenceCleanedDoi
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


## sampleWD
write.table(
  x = sampleWD,
  file = pathDataInterimTablesAnalysedSampleGoldONPDB,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
