# title: "integrating bio chemo ref"

# loading
##paths
source("paths.R")

## libraries
library(tidyverse)
library(dplyr)
library(readr)

## functions
source("functions/analysis.R")
source("functions/helpers.R")

## files
### original table
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### organism
organismTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedOrganismFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(!is.na(organismCleaned)) %>% # this step is new and should avoid useless too big files
  select(
    -organismCleaned,
    -organismDbTaxoQuality,
    -organism_1_kingdom,
    -organism_2_phylum,
    -organism_3_class,
    -organism_4_order,
    -organism_5_family,
    -organism_6_genus,
    -organism_7_species,
    -organism_8_variety
  ) %>%
  filter(!is.na(organismLowestTaxon)) # some DBs propose a taxon ID but with no name

### structure
#### translated
translatedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(!is.na(structureTranslated)) %>% # this step is new and should avoid useless too big files
  # mutate(
  #   structureOriginal_inchi = structureOriginalInchi, # remove this once reran
  #   structureOriginal_smiles = structureOriginalSmiles,
  #   structureOriginal_nominal = structureOriginalNominal,
  #   structureTranslated_smiles = structureTranslatedSmiles,
  #   structureTranslated_nominal = structureTranslatedNominal
  # ) %>%
  select(-structureTranslated_smiles,
         -structureTranslated_nominal,
         -nameCleaned) %>%
  distinct(
    structureOriginal_inchi,
    structureOriginal_smiles,
    structureOriginal_nominal,
    .keep_all = TRUE
  )

#### cleaned
cleanedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(!is.na(structureTranslated)) %>% # this step is new and should avoid useless too big files
  select(-validatorLog,
         -formulaSanitized,
         -exactmassSanitized,
         -xlogpSanitized)

### reference table
referenceTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  filter(
    !is.na(referenceOriginal_external) |
      !is.na(referenceOriginal_isbn) |
      !is.na(referenceCleanedTitle) |
      !is.na(referenceCleanedJournal) |
      !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmid) |
      !is.na(referenceCleanedPmcid)
  ) %>% # this step is new and should avoid useless too big files
  select(
    -referenceTranslatedDoi,
    -referenceTranslatedJournal,
    -referenceTranslatedTitle,
    -referenceTranslatedDate,
    -referenceTranslatedAuthor,
    -referenceTranslationScoreCrossref,
    -referenceTranslationScoreDistance,
    -referenceCleanedAuthor,
    -referenceCleanedDate
  )

# integrating
## structure
translatedStructureIntegratedTable <-
  left_join(originalTable, translatedStructureTable)

cleanedStructureIntegratedTable <-
  left_join(translatedStructureIntegratedTable, cleanedStructureTable)

## organism
organismStructureIntegratedTable <-
  left_join(cleanedStructureIntegratedTable, organismTable)

## reference
referenceOrganismStructureIntegratedTable <-
  left_join(organismStructureIntegratedTable, referenceTable)

# selecting minimal columns
inhouseDb <- referenceOrganismStructureIntegratedTable %>%
  mutate(
    referenceCleanedTranslationScoreCrossref = as.numeric(referenceCleanedTranslationScoreCrossref)
  )

# export
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesCurated),
  dir.create(pathDataInterimTablesCurated),
  FALSE
)

## table
write.table(
  x = inhouseDb,
  file = gzfile(
    description = pathDataInterimTablesCuratedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
