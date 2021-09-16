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
library(stringr)

log_debug("... functions")
source("r/capitalize.R")
source("r/name2inchi_cactus.R")
source("r/name2inchi_cts.R")
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
    "-ostdinchi",
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
  select(inchiNominal_opsin = `...1`)

dataInterim <- bind_cols(dataPreparedNamesDistinct, dataOpsin)

dataInterim <- left_join(dataPreparedNames, dataInterim)

log_debug("exporting interim ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal_opsin)
write_delim(
  x = dataInterim,
  file = pathDataInterimTablesTranslatedStructureNominal_opsin
)

dataForCTS <- dataInterim %>%
  filter(is.na(inchiNominal_opsin)) %>%
  distinct(nameCleaned, .keep_all = TRUE) %>%
  select(structureOriginal_nominal, nameCleaned)

log_debug("translating structures with CTS (slow but better results)")
if (nrow(dataForCTS) == 0) {
  dataForCTS[1, "nameCleaned"] <- NA
}

dataTranslatedNominal_cts <- dataForCTS %>%
  select(-structureOriginal_nominal) %>%
  mutate(inchiNominal_cts = invisible(
    pbmclapply(
      FUN = name2inchi_cts,
      X = seq_len(nrow(dataForCTS)),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = 2,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )) %>%
  mutate(inchiNominal_cts = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = inchiNominal_cts
  )) %>%
  mutate(inchiNominal_cts = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = inchiNominal_cts
  )) %>%
  mutate(inchiNominal_cts = y_as_na(
    x = inchiNominal_cts,
    y = "NA"
  ))

dataInterim_2 <- left_join(
  dataInterim,
  dataTranslatedNominal_cts
)

log_debug("exporting interim ...")
log_debug(
  pathDataInterimTablesTranslatedStructureNominal_cts
)
write_delim(
  x = dataInterim_2,
  file = pathDataInterimTablesTranslatedStructureNominal_cts
)

# dataInterim_2 <-
#   read_delim(file = pathDataInterimTablesTranslatedStructureNominal_cts)

dataInterim_2 <- dataInterim_2 %>%
  mutate(nameCleaned_capitalized = capitalize(nameCleaned))

dataForCTS_2 <- dataInterim_2 %>%
  filter(is.na(inchiNominal_cts)) %>%
  filter(nameCleaned_capitalized != nameCleaned)

log_debug("translating structures with CTS again (capitalized)")
if (nrow(dataForCTS_2) == 0) {
  dataForCTS_2[1, "nameCleaned_capitalized"] <- NA
}

dataTranslatedNominal_cts_2 <- dataForCTS_2 %>%
  select(-structureOriginal_nominal) %>%
  mutate(inchiNominal_cts_2 = invisible(
    pbmclapply(
      FUN = name2inchi_cts_capitalized,
      X = seq_len(nrow(dataForCTS_2)),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = 2,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )) %>%
  mutate(inchiNominal_cts_2 = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = inchiNominal_cts_2
  )) %>%
  mutate(inchiNominal_cts_2 = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = inchiNominal_cts_2
  )) %>%
  mutate(inchiNominal_cts_2 = y_as_na(
    x = inchiNominal_cts_2,
    y = "NA"
  ))

dataInterim_3 <- left_join(
  dataInterim_2,
  dataTranslatedNominal_cts_2
)

log_debug("exporting interim ...")
log_debug(
  pathDataInterimTablesTranslatedStructureNominal_cts_2
)
write_delim(
  x = dataInterim_3,
  file = pathDataInterimTablesTranslatedStructureNominal_cts_2
)

# dataInterim_3 <-
#   read_delim(file = pathDataInterimTablesTranslatedStructureNominal_cts_2)

## cactus is the lowest quality but allows retrieving important structures also
## some incorrect spotted...
### see https://cactus.nci.nih.gov/chemical/structure?identifier=terpinen-4-ol&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure?identifier=PONGAMOL&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure/Combretastatin%20b-2%20/smiles

dataForCactus <- dataInterim_3 %>%
  filter(is.na(inchiNominal_opsin)) %>%
  filter(is.na(inchiNominal_cts)) %>%
  filter(is.na(inchiNominal_cts_2))

if (nrow(dataForCactus) == 0) {
  dataForCactus[1, "nameCleaned"] <- NA
}

log_debug("... with cactus (fast)")
dataTranslatedNominal_cactus <- dataForCactus %>%
  select(-structureOriginal_nominal) %>%
  mutate(inchiNominal_cactus = invisible(
    pbmclapply(
      FUN = name2inchi_cactus,
      X = seq_len(nrow(dataForCactus)),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = 2,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )) %>%
  mutate(inchiNominal_cactus = as.character(inchiNominal_cactus)) %>%
  mutate(inchiNominal_cactus = y_as_na(
    x = inchiNominal_cactus,
    y = "character(0)"
  )) %>%
  mutate(inchiNominal_cactus = y_as_na(
    x = inchiNominal_cactus,
    y = "NA"
  )) %>%
  mutate(inchiNominal_cactus = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = inchiNominal_cactus
  )) %>%
  mutate(inchiNominal_cactus = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = inchiNominal_cactus
  ))

dataTranslated <- left_join(
  dataInterim_3,
  dataTranslatedNominal_cactus
) %>%
  mutate(inchiNominal_cts = ifelse(
    test = grepl(
      pattern = "^InChI=.*",
      x = inchiNominal_cts
    ),
    yes = inchiNominal_cts,
    no = NA
  )) %>%
  mutate(inchiNominal_cts_2 = ifelse(
    test = grepl(
      pattern = "^InChI=.*",
      x = inchiNominal_cts_2
    ),
    yes = inchiNominal_cts_2,
    no = NA
  )) %>%
  mutate(structureTranslated_nominal = ifelse(
    test = !is.na(inchiNominal_opsin),
    yes = inchiNominal_opsin,
    no = ifelse(
      test = !is.na(inchiNominal_cts),
      yes = inchiNominal_cts,
      no = ifelse(
        test = !is.na(inchiNominal_cts_2),
        yes = inchiNominal_cts_2,
        no = inchiNominal_cactus
      )
    )
  ))

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal)
write_delim(
  x = dataTranslated,
  file = pathDataInterimTablesTranslatedStructureNominal
)

# dataTranslated <-
#   read_delim(file = pathDataInterimTablesTranslatedStructureNominal)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
