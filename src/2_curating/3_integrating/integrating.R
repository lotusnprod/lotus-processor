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
    -smilesSanitized,
    -inchiSanitized,
    -formulaSanitized,
    -exactmassSanitized,
    -xlogpSanitized,
    -organismCleaned,
    -organismDbTaxoQuality,
    -organism_1_kingdom,
    -organism_2_phylum,
    -organism_3_class,
    -organism_4_order,
    -organism_5_family,
    -organism_6_genus,
    -organism_7_species,
    -organism_8_variety,
    -referenceTranslatedDoi,
    -referenceTranslatedJournal,
    -referenceTranslatedTitle,
    -referenceTranslatedDate,
    -referenceTranslatedAuthor,
    -referenceTranslationScore,
    -referenceCleanedAuthor,
    -referenceCleanedDate
  )

fullDbFiltered <- fullDb %>%
  filter(
    referenceCleanedTranslationScore >= 70 &
      referenceCleanedTranslationScore <= 150 &
      !is.na(organismCurated)
  ) %>%
  distinct(
    inchikeySanitized,
    organismCurated,
    referenceOriginalExternal,
    referenceCleanedDoi,
    .keep_all = TRUE
  )

fullDbFilteredDnp <- fullDb %>%
  filter(database == "dnp_1")

fullDbFilteredNoDnp <- fullDb %>%
  filter(database != "dnp_1")

fullDBDnpTop <- rbind(fullDbFilteredDnp, fullDbFilteredNoDnp)

fullDbFilteredOutsideDnp <- fullDBDnpTop %>%
  filter(!is.na(organismCurated)) %>%
  distinct(inchikeySanitized, organismCurated, .keep_all = TRUE) %>%
  filter(database != "dnp_1") %>%
  filter(referenceCleanedTranslationScore >= 70 &
           referenceCleanedTranslationScore <= 150)

stats <- fullDbFilteredOutsideDnp %>%
  group_by(database) %>%
  count()

## table
write.table(
  x = fullDb,
  file = gzfile(
    description = pathDataInterimTablesCleanedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
