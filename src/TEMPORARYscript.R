# temporary file to get some numbers to see if we are on the right way!

# loading
## functions
source("functions/analysis.R")
source("functions/helpers.R")

## libraries
library(dplyr)
library(readr)
library(tidyverse)

## files
### original table
originalTable <- read_delim(
  file = gzfile(description = "../data/interim/tables/0_original/originalTable.tsv.zip"),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### structure table
#### loading multiple old files, will be optimized later on
smilesStructureTable <- read_delim(
  file = gzfile(description = "../data/interim/tables/1_translated/translatedStructureSmiles.tsv.zip"),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

nominalStructureTable <- read_delim(
  file = gzfile(description = "../data/interim/tables/1_translated/translatedStructureNominal.tsv.zip"),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#### this is an old version of PM, might be that we got some very small differences. Do not check path, names and so on
structureTable <- read_delim(
  file = "../data/interim/tables/2_cleaned/opennpdbSanitizedClassyfied.tsv",
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### organism table
organismTable <- read_delim(
  file = gzfile(description = "../data/interim/tables/3_curated/curatedOrganism.tsv.zip"),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### reference table
referenceTable <- read_delim(
  file = gzfile(description = "../data/interim/tables/2_cleaned/cleanedReference.tsv.zip"),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting adequate minimal columns
## structure
### modifying some column names to have them closer to what they should look like in the end
structureTable <- structureTable %>%
  select(
    structureTranslated,
    smilesSanitized = smiles_sanitized,
    inchiSanitized = inchi_sanitized,
    inchiKeySanitized = inchikeySanitized,
    inchiKey2DSanitized = shortikSanitized,
    formulaSanitized,
    exactMassSanitized = exactmassSanitized,
    xlogpSanitized,
    structure_1_kingdom = kingdom.name,
    structure_2_superclass = superclass.name,
    structure_3_class = class.name,
    structure_4_subclass = subclass.name,
    structure_5_directParent = direct_parent.name
  )

## organism
organismTable <- organismTable %>%
  select(...)

## organism
referenceTable <- referenceTable %>%
  select(...)

# integrating
## structure
### some additional manipulations to be done because of merging of old data but nothing terrible
structureIntegratedTable <-
  left_join(originalTable, smilesStructureTable)

structureIntegratedTable <-
  left_join(structureIntegratedTable, nominalStructureTable) %>%
  select(
    structureOriginalInchi,
    structureTranslatedSmiles,
    structureTranslatedNominal,
    everything()
  )

structureIntegratedTable$structureTranslated <-
  apply(structureIntegratedTable[1:3],
        1,
        function(x)
          tail(na.omit(x), 1))

structureIntegratedTable$structureTranslated <-
  as.character(structureIntegratedTable$structureTranslated)

structureIntegratedTable$structureTranslated <-
  y_as_na(x = structureIntegratedTable$structureTranslated,
          y = "character(0)")

structureIntegratedTable$structureTranslated <-
  y_as_na(x = structureIntegratedTable$structureTranslated,
          y = "NA")

structureIntegratedTable <-
  left_join(structureIntegratedTable, structureTable)

## organism
organismStructureIntegratedTable <-
  left_join(structureIntegratedTable, organismTable)

## reference
referenceOrganismStructureIntegratedTable <-
  left_join(organismStructureIntegratedTable, referenceTable)

# selecting minimal columns
fullDb <- referenceOrganismStructureIntegratedTable %>%
  mutate(cleanedTranslationScore = as.numeric(cleanedTranslationScore)) %>% 
  select(database,
         organismOriginal,
         structureTranslated,
         referenceOriginal,
         inchiSanitized,
         inchiKeySanitized,
         inchiKey2DSanitized,
         organismCurated,
         cleanedDoi,
         cleanedTranslationScore
         )

fullDbFiltered <- fullDb %>% 
  filter(cleanedTranslationScore >= 90 &
           cleanedTranslationScore <= 110 &
           !is.na(organismCurated)) %>% 
  distinct(inchiSanitized, organismCurated, cleanedDoi, .keep_all = TRUE)
