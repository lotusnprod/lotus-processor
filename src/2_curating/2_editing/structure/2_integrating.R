cat("This script integrates all chemical translations. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
source("r/vroom_safe.R")

cat("loading files ... \n")
cat("... whole chemicals list \n")
originalTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalStructureFull)

cat("... chemical names list \n")
nominalStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedStructureNominal)

cat("... SMILES list \n")
smilesStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedStructureSmiles)

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

if (nrow(translatedStructureTable) == 0) {
  translatedStructureTable[1, c(
    "structureType",
    "structureValue",
    "structureTranslated"
  )] <- NA
}

cat("outputing unique structures \n")
translatedStructureTableUnique <- translatedStructureTable %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

if (nrow(translatedStructureTableUnique) == 0) {
  translatedStructureTableUnique[1, "structureTranslated"] <- NA
}

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedStructureFinal, "\n")
vroom_write(
  x = translatedStructureTable,
  path = pathDataInterimTablesTranslatedStructureFinal
)

cat(pathDataInterimTablesTranslatedStructureUnique, "\n")
vroom_write(
  x = translatedStructureTableUnique,
  path = pathDataInterimTablesTranslatedStructureUnique
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
