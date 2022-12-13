source("r/log_debug.R")
log_debug("This script performs chemical name translation.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(readr)
library(splitstackshape)
library(stringr)

log_debug("... functions")
source("r/capitalize.R")
source("r/name2structure_cactus.R")
source("r/name2structure_pubchem.R")
source("r/preparing_name.R")
source("r/progressr.R")
source("r/y_as_na.R")

log_debug("loading chemical names lists")
dataOriginal <-
  readr::read_delim(
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

dataPreparedNamesDistinct <- dataPreparedNames |>
  dplyr::distinct(nameCleaned)

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimTablesTranslatedStructure)

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
  readr::read_delim(
    file = pathDataInterimTablesTranslatedStructureOpsin,
    delim = "\t",
    skip_empty_rows = FALSE,
    escape_double = TRUE,
    trim_ws = TRUE
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::select(smilesNominal_opsin = `...1`)

dataInterim <-
  dplyr::bind_cols(dataPreparedNamesDistinct, dataOpsin)

dataInterim <- dplyr::left_join(dataPreparedNames, dataInterim) |>
  dplyr::distinct()

log_debug("exporting interim ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal_opsin)
readr::write_delim(
  x = dataInterim,
  file = pathDataInterimTablesTranslatedStructureNominal_opsin,
  delim = "\t",
  na = ""
)

dataForPubchem <- dataInterim |>
  dplyr::filter(is.na(smilesNominal_opsin)) |>
  dplyr::distinct(nameCleaned)

if (nrow(dataForPubchem) == 0) {
  dataForPubchem[1, "nameCleaned"] <- NA
}

log_debug("translating structures with pubchem")
dataTranslatedNominal_pubchem <- dataForPubchem |>
  dplyr::mutate(smilesNominal_pubchem = name2smiles_pubchem(xs = seq_len(
    nrow(dataForPubchem)
  ) |>
    progressr::with_progress())) |>
  splitstackshape::cSplit("smilesNominal_pubchem",
    sep = "\n",
    direction = "long"
  ) |>
  dplyr::mutate(smilesNominal_pubchem = y_as_na(
    x = trimws(smilesNominal_pubchem),
    y = "NA"
  )) |>
  dplyr::mutate(
    smilesNominal_pubchem = gsub(
      pattern = "^http.*",
      replacement = NA,
      x = smilesNominal_pubchem
    )
  )

dataInterim_1 <- dplyr::left_join(
  dataInterim,
  dataTranslatedNominal_pubchem
) |>
  dplyr::distinct()

log_debug("exporting interim ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal_pubchem)
readr::write_delim(
  x = dataInterim_1,
  file = pathDataInterimTablesTranslatedStructureNominal_pubchem,
  delim = "\t",
  na = ""
)

# dataInterim_1 <- readr::read_delim(
#   file = pathDataInterimTablesTranslatedStructureNominal_pubchem,
#   delim = "\t"
# )

## cactus is the lowest quality but allows retrieving important structures also
## some incorrect spotted...
### see https://cactus.nci.nih.gov/chemical/structure?identifier=terpinen-4-ol&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure?identifier=PONGAMOL&representation=smiles
### or https://cactus.nci.nih.gov/chemical/structure/Combretastatin%20b-2%20/smiles

dataForCactus <- dataInterim_1 |>
  dplyr::filter(is.na(smilesNominal_opsin)) |>
  dplyr::filter(is.na(smilesNominal_pubchem)) |> 
  dplyr::select(-structureOriginal_nominal) |>
  dplyr::distinct()

if (nrow(dataForCactus) == 0) {
  dataForCactus[1, "nameCleaned"] <- NA
}

log_debug("... with cactus (fast)")
dataTranslatedNominal_cactus <- dataForCactus |>
  dplyr::mutate(smilesNominal_cactus = name2smiles_cactus(xs = seq_len(
    nrow(dataForCactus)
  ) |>
    progressr::with_progress())) |>
  dplyr::mutate(smilesNominal_cactus = as.character(smilesNominal_cactus)) |>
  dplyr::mutate(smilesNominal_cactus = y_as_na(
    x = smilesNominal_cactus,
    y = "character(0)"
  )) |>
  dplyr::mutate(smilesNominal_cactus = y_as_na(
    x = smilesNominal_cactus,
    y = "NA"
  )) |>
  dplyr::mutate(smilesNominal_cactus = gsub(
    pattern = "^http.*",
    replacement = NA,
    x = smilesNominal_cactus
  )) |>
  dplyr::mutate(smilesNominal_cactus = gsub(
    pattern = "^NCI.*",
    replacement = NA,
    x = smilesNominal_cactus
  ))

dataTranslated <- dplyr::left_join(
  dataInterim_1,
  dataTranslatedNominal_cactus
) |>
  dplyr::mutate(structureTranslated_nominal = ifelse(
    test = !is.na(smilesNominal_opsin),
    yes = smilesNominal_opsin,
    no = ifelse(
      test = !is.na(smilesNominal_pubchem),
      yes = smilesNominal_pubchem,
      no = smilesNominal_cactus
    )
  )) |>
  dplyr::distinct()

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedStructureNominal)
readr::write_delim(
  x = dataTranslated,
  file = pathDataInterimTablesTranslatedStructureNominal,
  delim = "\t",
  na = ""
)

# dataTranslated <-
#   readr::read_delim(file = pathDataInterimTablesTranslatedStructureNominal,
#              delim = "\t")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
