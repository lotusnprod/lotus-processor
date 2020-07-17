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
)

### structure
#### translated
translatedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  distinct(
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    .keep_all = TRUE
  )

#### cleaned
cleanedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### reference table
referenceTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting adequate minimal columns
## structure
# structureTable <- structureTable %>%
#   select(...)

# ## organism
# organismTable <- organismTable %>%
#   select(...)

# ## organism
# referenceTable <- referenceTable %>%
#   select(...)

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
fullDb <- referenceOrganismStructureIntegratedTable %>%
  mutate(referenceCleanedTranslationScore = as.numeric(referenceCleanedTranslationScore)) %>%
  select(
    -structureTranslatedSmiles,
    -structureTranslatedNominal,
    -nameCleaned,
    -structureTranslated,
    -validatorLog,
    -formulaSanitized,
    -exactmassSanitized,
    -xlogpSanitized,
    -organismCleaned,
    -organismDbTaxoQuality,
    -referenceTranslatedDoi,
    -referenceTranslatedJournal,
    -referenceTranslatedTitle,
    -referenceTranslatedDate,
    -referenceTranslatedAuthor,
    -referenceTranslationScore,
    -referenceCleanedAuthor,
    -referenceCleanedDate
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
  x = fullDb,
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
