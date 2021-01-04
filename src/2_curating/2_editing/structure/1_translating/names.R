cat("This script performs chemical name translation. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(pbmcapply)

cat("... functions \n")
source("r/y_as_na.R")
source("r/preparing_name.R")
source("r/name2inchi_cactus.R")
source("r/name2inchi_cts.R")
source("r/vroom_safe.R")

cat("loading chemical names lists \n")
dataOriginal <-
  vroom_read_safe(path = pathDataInterimTablesOriginalStructureNominal)

if (nrow(dataOriginal) == 0) {
  dataOriginal[1, "structureOriginal_nominal"] <- NA
}

# to avoid errors if dataframe not empty at the begining but filled with NA
colnames(dataOriginal)[1] <- "structureOriginal_nominal"

cat("preparing names \n")
dataPreparedNames <- preparing_name(x = dataOriginal)

dataPreparedNamesDistinct <- dataPreparedNames %>%
  distinct(nameCleaned)

cat("ensuring directories exist \n")
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

cat("exporting prepared names ... \n")
cat(pathDataInterimTablesTranslatedStructureNominal, "\n")
write.table(
  x = dataPreparedNamesDistinct,
  file = pathDataInterimTablesTranslatedStructurePrepared_1,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat("translating names ... \n")
cat("... with OPSIN (very fast but limited) \n")

system(
  command = paste(
    "java -jar",
    pathBinOpsin,
    "-ostdinchi",
    pathDataInterimTablesTranslatedStructurePrepared_1,
    pathDataInterimTablesTranslatedStructureOpsin
  )
)

cat("loading opsin results \n")
dataOpsin <-
  read_delim(
    file = pathDataInterimTablesTranslatedStructureOpsin,
    delim = "\t",
    skip_empty_rows = FALSE,
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  mutate_all(as.character) %>%
  select(inchiNominal_opsin = X1)

dataInterim <- bind_cols(dataPreparedNamesDistinct, dataOpsin)

dataInterim <- left_join(dataPreparedNames, dataInterim)

cat("exporting interim ... \n")
cat(pathDataInterimTablesTranslatedStructureNominal_opsin, "\n")
vroom_write_safe(
  x = dataInterim,
  path = pathDataInterimTablesTranslatedStructureNominal_opsin
)

dataForCactus <- dataInterim %>%
  filter(is.na(inchiNominal_opsin)) %>%
  distinct(nameCleaned, .keep_all = TRUE) %>%
  select(structureOriginal_nominal, nameCleaned)

## some incorrect spotted...
### see https://cactus.nci.nih.gov/chemical/structure?identifier=terpinen-4-ol&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure?identifier=PONGAMOL&representation=smiles

cat("... with cactus (fast) \n")
dataTranslatedNominal_cactus <- dataForCactus %>%
  mutate(inchiNominal_cactus = invisible(
    pbmclapply(
      FUN = name2inchi_cactus,
      X = 1:nrow(dataForCactus),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
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

dataInterim_2 <- left_join(
  dataInterim,
  dataTranslatedNominal_cactus
)

cat("exporting interim ... \n")
cat(
  pathDataInterimTablesTranslatedStructureNominal_cactus,
  "\n"
)
vroom_write_safe(
  x = dataInterim_2,
  path = pathDataInterimTablesTranslatedStructureNominal_cactus
)

dataForCTS <- dataInterim_2 %>%
  filter(is.na(inchiNominal_opsin) & is.na(inchiNominal_cactus)) %>%
  distinct(nameCleaned, .keep_all = TRUE) %>%
  select(structureOriginal_nominal, nameCleaned)

cat("translating structures with CTS (slow but more results) \n")
if (nrow(dataForCTS) == 0) {
  dataForCTS[1, "nameCleaned"] <- NA
}

dataTranslatedNominal_cts <- dataForCTS %>%
  mutate(inchiNominal_cts = invisible(
    pbmclapply(
      FUN = name2inchi_cts,
      X = 1:nrow(dataForCTS),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
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

dataTranslated <- left_join(
  dataInterim_2,
  dataTranslatedNominal_cts
) %>%
  mutate(inchiNominal_cts = ifelse(
    test = grepl(
      pattern = "^InChI=.*",
      x = inchiNominal_cts
    ),
    yes = inchiNominal_cts,
    no = NA
  )) %>%
  mutate(structureTranslated_nominal = ifelse(
    test = !is.na(inchiNominal_opsin),
    yes = inchiNominal_opsin,
    no = ifelse(
      test = !is.na(inchiNominal_cactus),
      yes = inchiNominal_cactus,
      no = inchiNominal_cts
    )
  ))

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedStructureNominal, "\n")
vroom_write_safe(
  x = dataTranslated,
  path = pathDataInterimTablesTranslatedStructureNominal
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
