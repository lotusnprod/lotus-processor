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
  read_delim(
    file = pathDataInterimTablesOriginalStructureFull,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... chemical names list")
nominalStructureTable <-
  read_delim(
    file = pathDataInterimTablesTranslatedStructureNominal,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... SMILES list")
inchiStructureTable <-
  read_delim(
    file = pathDataInterimTablesTranslatedStructureInchi,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("joining")
translatedStructureTable <-
  left_join(
    originalTable,
    inchiStructureTable,
    by = c("structureValue" = "structureOriginal_inchi")
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
      test = structureType == "smiles",
      yes = structureValue,
      no = ifelse(
        test = structureType == "inchi",
        yes = structureTranslated_inchi,
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

log_debug("filtering out bad SMILES for next step")
## This list will be outsourced later on
bad_smiles <- c(
  "[OH-].[OH-].[OH-].[OH-].[OH-].[OH-].[V+6]", ## rdkit fails
  "Not Available",
  "Cc1ccc2c(c1)-n1-c(=O)/c=c\\c(=O)-n-2-c2cc(C)ccc2-1", ## in reality 'Cc1ccc2c(c1)-n1-c(=O)/c=c\c(=O)-n-2-c2cc(C)ccc2-1'
  "CCCCCc1cccc([O-])c1C1=NC(C2SCC(C([O-])C(C)(C)C3=NC(C)(C(=O)[O-])CS3)N2C)CS1.[Zn+3]",
  "CC(C)C([NH+]=CO)C(=O)OC(C)C1C/C=C/CCC(O)C(C)(C)C2CCC(C)C3(O2)O[B-]24OC(C(=O)O1)C1(OC(CCC1C)C(C)(C)C(O)CCCC1CC(OC(=O)C3O2)C(C)O1)O4",
  "CC(O)=[NH+]C(C(=O)OC(C)C1C/C=C/CCC(O)C(C)(C)C2CCC(C)C3(O2)O[B-]24OC(C(=O)O1)C1(OC(CCC1C)C(C)(C)C(O)CCCC1CC(OC(=O)C3O2)C(C)O1)O4)C(C)C"
  # "CC(C)c1c(O)cc(CCc2ccccc2)cc1OC(=N)O", # primary carbamate fail
  # "COc1c(N=C(O)c2ccc(NC(=O)c3ccc(NC(=O)C(NC(=O)c4ccc(NC(=O)/C(C)=C/c5ccc(OC(=N)O)cc5)cc4)C(C#N)OC)cc3)c(OC)c2O)ccc(C(=O)O)c1O", # primary carbamate fail
  # "COc1c(N=C(O)c2ccc(NC(=O)c3ccc(NC(=O)C(NC(=O)c4ccc(NC(=O)/C(C)=C/c5ccc(OC(=N)O)cc5)cc4)C(OC)C(=N)O)cc3)c(OC)c2O)ccc(C(=O)O)c1O" # primary carbamate fail
)

translatedStructureTable <- translatedStructureTable %>%
  filter(!structureTranslated %in% bad_smiles) %>%
  mutate(di = gsub(
    pattern = "\\..*",
    replacement = "",
    x = structureTranslated
  )) %>% filter(
    structureTranslated == di | gsub(
      pattern = ".*\\.",
      replacement = "",
      x = structureTranslated
    ) == di
  ) %>%
  select(structureType, 
         structureValue, 
         structureTranslated = di) %>%
  distinct()

log_debug("outputing unique structures")
translatedStructureTableUnique <- translatedStructureTable %>%
  distinct(structureTranslated)

if (nrow(translatedStructureTableUnique) == 0) {
  translatedStructureTableUnique[1, "structureTranslated"] <- NA
}

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedStructureFinal)
write_delim(
  x = translatedStructureTable,
  file = pathDataInterimTablesTranslatedStructureFinal,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimTablesTranslatedStructureUnique)
write_delim(
  x = translatedStructureTableUnique,
  file = pathDataInterimTablesTranslatedStructureUnique,
  delim = "\t",
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
