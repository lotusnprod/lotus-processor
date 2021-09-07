source("r/log_debug.R")
log_debug("This script integrates all chemical translations.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("loading files ...")
log_debug("... whole chemicals list")
originalTable <-
  read_delim(file = pathDataInterimTablesOriginalStructureFull,
             delim = "\t")

log_debug("... chemical names list")
nominalStructureTable <-
  read_delim(file = pathDataInterimTablesTranslatedStructureNominal)

log_debug("... SMILES list")
smilesStructureTable <-
  read_delim(file = pathDataInterimTablesTranslatedStructureSmiles)

log_debug("joining")
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

log_debug("outputing unique structures")
translatedStructureTableUnique <- translatedStructureTable %>%
  filter(!is.na(structureTranslated)) %>%
  distinct(structureTranslated)

if (nrow(translatedStructureTableUnique) == 0) {
  translatedStructureTableUnique[1, "structureTranslated"] <- NA
}

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedStructureFinal)
write_delim(
  x = translatedStructureTable,
  file = pathDataInterimTablesTranslatedStructureFinal
)

log_debug(pathDataInterimTablesTranslatedStructureUnique)
write_delim(
  x = translatedStructureTableUnique,
  file = pathDataInterimTablesTranslatedStructureUnique
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
