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
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>% 
  select(structureOriginalInchi,
         structureOriginalSmiles,
         structureOriginalNominal)

### structure table
#### loading multiple old files, will be optimized later on
smilesStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedSmiles),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

nominalStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedNominal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# joining
translatedStructureTable <-
  left_join(originalTable, smilesStructureTable)
translatedStructureTable <-
  left_join(translatedStructureTable, nominalStructureTable) %>%
  select(
    structureOriginalInchi,
    structureTranslatedSmiles,
    structureTranslatedNominal,
    everything()
  )

translatedStructureTable$structureTranslated <-
  apply(translatedStructureTable[1:3],
        1,
        function(x)
          tail(na.omit(x), 1))

translatedStructureTable$structureTranslated <-
  as.character(translatedStructureTable$structureTranslated)

translatedStructureTable$structureTranslated <-
  y_as_na(x = translatedStructureTable$structureTranslated,
          y = "character(0)")

translatedStructureTable$structureTranslated <-
  y_as_na(x = translatedStructureTable$structureTranslated,
          y = "NA")

# unique structures
translatedStructureTableUnique <- translatedStructureTable %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

# exporting
write.table(
  x = translatedStructureTable,
  file = gzfile(
    description = pathDataInterimTablesTranslatedStructure,
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
