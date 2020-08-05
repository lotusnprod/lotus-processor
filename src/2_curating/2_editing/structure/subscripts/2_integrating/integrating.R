# title: "integrating chemo"

# loading
## paths
source("paths.R")

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
  file = gzfile(description = pathDataInterimTablesOriginalStructureFull),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

### structure table
#### loading multiple old files, will be optimized later on
nominalStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureNominal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

smilesStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureSmiles),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# joining
translatedStructureTable <-
  left_join(
    originalTable,
    smilesStructureTable,
    by = c("structureValue" = "structureOriginal_smiles")
  )

translatedStructureTable <-
  left_join(
    translatedStructureTable,
    nominalStructureTable,
    by = c("structureValue" = "structureOriginal_nominal")
  )

translatedStructureTable <- translatedStructureTable %>%
  mutate(
    structureTranslated = ifelse(
      test = structureType == "inchi",
      yes = structureValue,
      no = ifelse(
        test = structureType == "smiles",
        yes = structureTranslated_smiles,
        no = structureTranslated_nominal
      )
    )
  ) %>%
  distinct(structureType, structureValue, structureTranslated) %>%
  filter(!is.na(structureTranslated))

# unique structures
translatedStructureTableUnique <- translatedStructureTable %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

# exporting
write.table(
  x = translatedStructureTable,
  file = gzfile(
    description = pathDataInterimTablesTranslatedStructureFinal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

write.table(
  x = translatedStructureTableUnique,
  file = gzfile(
    description = pathDataInterimTablesTranslatedStructureUnique,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
