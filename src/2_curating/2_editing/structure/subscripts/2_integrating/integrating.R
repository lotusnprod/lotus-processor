cat("This script integrates all chemical translations. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/analysis.R")
source("functions/helpers.R")

cat("loading files ... \n")
cat("... whole chemicals list \n")
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalStructureFull),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("... chemical names list \n")
nominalStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureNominal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("... SMILES list \n")
smilesStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureSmiles),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("joining \n")
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

if (nrow(translatedStructureTable) == 0)
  translatedStructureTable[1,] <- NA

cat("outputing unique structures \n")
translatedStructureTableUnique <- translatedStructureTable %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

if (nrow(translatedStructureTableUnique) == 0)
  translatedStructureTableUnique[1,] <- NA

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedStructureFinal, "\n")
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

cat(pathDataInterimTablesTranslatedStructureUnique, "\n")
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

end <- Sys.time()

cat("Script finished in", end - start , "seconds \n")
