source("r/log_debug.R")
log_debug("This script performs chemical name translation.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(readr)
library(splitstackshape)
library(stringr)

log_debug("... functions")
source("r/capitalize.R")
source("r/name2structure_cactus.R")
source("r/name2structure_pubchem.R")
source("r/parallel.R")
source("r/preparing_name.R")
source("r/y_as_na.R")

log_debug("loading chemical names lists")
dataOriginal <-
  read_delim(
    file = pathDataInterimTablesOriginalStructureNominal,
    delim = "\t"
  )

if (nrow(dataOriginal) == 0) {
  dataOriginal[1, "structureOriginal_nominal"] <- NA
}

## to avoid errors if dataframe not empty at the begining but filled with NA
colnames(dataOriginal)[1] <- "structureOriginal_nominal"

log_debug("preparing names")
dataPreparedNames <- preparing_name(x = dataOriginal)

dataPreparedNamesDistinct <- dataPreparedNames %>%
  distinct(nameCleaned)

log_debug("ensuring directories exist")
ifelse(
  test = !dir.exists(pathDataInterimTablesTranslated),
  yes = dir.create(pathDataInterimTablesTranslated),
  no = paste(pathDataInterimTablesTranslated, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesTranslatedStructure),
  yes = dir.create(pathDataInterimTablesTranslatedStructure),
  no = paste(pathDataInterimTablesTranslatedStructure, "exists")
)

log_debug("exporting prepared names ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal)

write.table(
  x = dataPreparedNamesDistinct,
  file = pathDataInterimTablesTranslatedStructurePrepared_1,
  na = "",
  row.names = FALSE,
  quote = FALSE,
  qmethod = "double",
  sep = "\t",
  fileEncoding = "UTF-8"
)

log_debug("translating names ...")
log_debug("... with OPSIN (very fast but limited)")

system(
  command = paste(
    "java -jar",
    pathBinOpsin,
    "-o smi",
    pathDataInterimTablesTranslatedStructurePrepared_1,
    pathDataInterimTablesTranslatedStructureOpsin
  )
)

log_debug("loading opsin results")
dataOpsin <-
  read_delim(
    file = pathDataInterimTablesTranslatedStructureOpsin,
    delim = "\t",
    skip_empty_rows = FALSE,
    escape_double = TRUE,
    trim_ws = TRUE
  ) %>%
  mutate_all(as.character) %>%
  select(smilesNominal_opsin = `...1`)

dataInterim <- bind_cols(dataPreparedNamesDistinct, dataOpsin)

dataInterim <- left_join(dataPreparedNames, dataInterim) %>%
  distinct()

log_debug("exporting interim ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal_opsin)
write_delim(
  x = dataInterim,
  file = pathDataInterimTablesTranslatedStructureNominal_opsin,
  delim = "\t",
  na = ""
)

dataForPubchem <- dataInterim %>%
  filter(is.na(smilesNominal_opsin)) %>%
  distinct(nameCleaned, .keep_all = TRUE) %>%
  select(structureOriginal_nominal, nameCleaned)

if (nrow(dataForPubchem) == 0) {
  dataForPubchem[1, "nameCleaned"] <- NA
}

log_debug("translating structures with pubchem")
dataTranslatedNominal_pubchem <- dataForPubchem %>%
  select(-structureOriginal_nominal) %>%
  mutate(smilesNominal_pubchem = invisible(
    pbmclapply(
      FUN = name2smiles_pubchem,
      X = seq_len(nrow(dataForPubchem)),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.cores = if (.Platform$OS.type == "unix") {
        2 ## limit to 5 calls per sec
      } else {
        1
      },
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
  )) %>%
  cSplit("smilesNominal_pubchem",
    sep = "\n",
    direction = "long"
  ) %>%
  mutate(smilesNominal_pubchem = y_as_na(
    x = trimws(smilesNominal_pubchem),
    y = "NA"
  )) %>%
  mutate(smilesNominal_pubchem = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = smilesNominal_pubchem
  ))

dataInterim_1 <- left_join(
  dataInterim,
  dataTranslatedNominal_pubchem
) %>%
  distinct()

log_debug("exporting interim ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal_pubchem)
write_delim(
  x = dataInterim_1,
  file = pathDataInterimTablesTranslatedStructureNominal_pubchem,
  delim = "\t",
  na = ""
)

# dataInterim_1 <- read_delim(
#   file = pathDataInterimTablesTranslatedStructureNominal_pubchem,
#   delim = "\t"
# )

## cactus is the lowest quality but allows retrieving important structures also
## some incorrect spotted...
### see https://cactus.nci.nih.gov/chemical/structure?identifier=terpinen-4-ol&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure?identifier=PONGAMOL&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure/Combretastatin%20b-2%20/smiles

dataForCactus <- dataInterim_1 %>%
  filter(is.na(smilesNominal_opsin)) %>%
  filter(is.na(smilesNominal_pubchem))

if (nrow(dataForCactus) == 0) {
  dataForCactus[1, "nameCleaned"] <- NA
}

log_debug("... with cactus (fast)")
dataTranslatedNominal_cactus <- dataForCactus %>%
  select(-structureOriginal_nominal) %>%
  mutate(smilesNominal_cactus = invisible(
    pbmclapply(
      FUN = name2smiles_cactus,
      X = seq_len(nrow(dataForCactus)),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.cores = numCores,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
  )) %>%
  mutate(smilesNominal_cactus = as.character(smilesNominal_cactus)) %>%
  mutate(smilesNominal_cactus = y_as_na(
    x = smilesNominal_cactus,
    y = "character(0)"
  )) %>%
  mutate(smilesNominal_cactus = y_as_na(
    x = smilesNominal_cactus,
    y = "NA"
  )) %>%
  mutate(smilesNominal_cactus = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = smilesNominal_cactus
  )) %>%
  mutate(smilesNominal_cactus = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = smilesNominal_cactus
  ))

dataTranslated <- left_join(
  dataInterim_1,
  dataTranslatedNominal_cactus
) %>%
  mutate(structureTranslated_nominal = ifelse(
    test = !is.na(smilesNominal_opsin),
    yes = smilesNominal_opsin,
    no = ifelse(
      test = !is.na(smilesNominal_pubchem),
      yes = smilesNominal_pubchem,
      no = smilesNominal_cactus
    )
  )) %>%
  distinct()

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal)
write_delim(
  x = dataTranslated,
  file = pathDataInterimTablesTranslatedStructureNominal,
  delim = "\t",
  na = ""
)

# dataTranslated <-
#   read_delim(file = pathDataInterimTablesTranslatedStructureNominal,
#              delim = "\t")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
