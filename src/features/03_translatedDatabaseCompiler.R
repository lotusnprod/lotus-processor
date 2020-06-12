# title: "Open NP DB (translated) compileR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
## original
dataOriginal <- read_delim(
  file = gzfile(pathOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## organism
dataTranslatedOrganism <- read_delim(
  file = gzfile(pathTranslatedOrganism),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## inchi
dataTranslatedStructureInchi <- read_delim(
  file = gzfile(pathOriginalStructureInchi),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## smiles
dataTranslatedStructureSmiles <- read_delim(
  file = gzfile(pathTranslatedStructureSmiles),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## nominal
dataTranslatedStructureNominal <- read_delim(
  file = gzfile(pathTranslatedStructureNominal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## references
dataTranslatedReference <- read_delim(
  file = gzfile(pathTranslatedReference),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# preparing data
## chemo
### inchi
dataTranslated <-
  left_join(dataOriginal,
            dataTranslatedStructureInchi)

### smiles
dataTranslated <-
  left_join(dataTranslated,
            dataTranslatedStructureSmiles)

### nominal
dataTranslated <-
  left_join(dataTranslated,
            dataTranslatedStructureNominal)

# arranging column order
dataTranslated <- dataTranslated %>%
  select(
    database,
    name,
    nameSanitized,
    organismOriginal,
    referenceOriginal,
    structureOriginalSmiles,
    structureOriginalNominal,
    structureOriginalInchi,
    structureTranslatedSmiles,
    structureTranslatedNominal,
  )

# single translated inchi column
dataTranslated$structureTranslated <-
  apply(dataTranslated[(ncol(dataTranslated) - 2):ncol(dataTranslated)],
        1,
        function(x)
          tail(na.omit(x), 1))

dataTranslated$structureTranslated <-
  as.character(dataTranslated$structureTranslated)

dataTranslated$structureTranslated <-
  y_as_na(x = dataTranslated$structureTranslated,
          y = "character(0)")

dataTranslated$structureTranslated <-
  y_as_na(x = dataTranslated$structureTranslated,
          y = "NA")

dataTranslated$structureTranslated <-
  y_as_na(x = dataTranslated$structureTranslated,
          y = "N/A")

## reference
dataTranslated <-
  left_join(dataTranslated,
            dataTranslatedReference)

## bio
dataTranslatedFilled <-
  left_join(dataTranslated,
            dataTranslatedOrganism)

# selecting
## structure table
dataTranslatedStructure <- dataTranslatedFilled %>%
  select(
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    structureTranslated
  ) %>%
  distinct(
    structureOriginalInchi,
    structureOriginalSmiles,
    structureOriginalNominal,
    .keep_all = TRUE
  )

## distinct structures for rdkit
dataTranslatedStructureRdkit <- dataTranslatedFilled %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

## distinct organisms for gnfinder
dataTranslatedOrganismGnfinder <- dataTranslatedFilled %>%
  filter(!is.na(organismTranslated)) %>%
  distinct(organismTranslated) %>%
  data.table()

# Chose next step for references
dataTranslatedReference <- dataTranslatedFilled %>%
  select(
    referenceOriginal,
    referenceSplit,
    translatedTitle,
    translatedJournal,
    translatedDate,
    translatedAuthor,
    translatedDoi,
    translationScore
  ) %>%
  distinct(
    referenceOriginal,
    referenceSplit,
    translatedTitle,
    translatedJournal,
    translatedDoi,
    .keep_all = TRUE
  )

# exporting
## structures for rdkit
write.table(
  x = dataTranslatedStructureRdkit,
  file = gzfile(
    description = pathTranslatedStructureDistinct,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## organisms for gnfinder
split_data_table(x = dataTranslatedOrganismGnfinder,
                 no_rows_per_frame = 10000,
                 path_to_store = pathTranslatedOrganismDistinct)

## references
write.table(
  x = dataTranslatedReference,
  file = gzfile(
    description = pathTranslatedReference,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

## table
write.table(
  x = dataTranslatedFilled,
  file = gzfile(
    description = pathTranslatedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)