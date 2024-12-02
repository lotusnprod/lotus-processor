source("r/log_debug.R")
log_debug("This script adds np-classifier classes to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(data.table)
library(dplyr)
library(future)
library(future.apply)
library(jsonlite)
library(progressr)
library(RCurl)
library(readr)
library(tidyr)

log_debug("... functions")
source("r/getClass.R")
source("r/progressr.R")
source("r/treat_npclassifier_json.R")

log_debug("loading smiles ...")
smiles <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedStructureStereoCounted,
    delim = "\t",
    col_types = cols(.default = "c")
  ) |>
  dplyr::distinct(structure_smiles_2D = smilesSanitizedFlat)
# smiles <-
#   readr::read_delim(
#     file = pathDataInterimDictionariesStructureMetadata,
#     delim = "\t",
#     col_types = cols(.default = "c")
#   ) |>
#   dplyr::distinct(structure_smiles_2D = structureCleaned_smiles2D)

log_debug("loading npClassifier taxonomy ...")
taxonomy <- jsonlite::fromJSON(txt = list.files(
  path = file.path(
    pathDataExternal,
    "taxonomySource/structure/npclassifier/"
  ),
  pattern = ".json",
  full.names = TRUE
))

if (file.exists(pathDataInterimDictionariesStructureDictionaryNpclassifierFile)) {
  old <-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile,
      delim = "\t",
      col_types = cols(.default = "c")
    ) |>
    dplyr::distinct() |>
    dplyr::mutate_all(as.character)

  new <- dplyr::anti_join(smiles, old)
} else {
  new <- smiles
}

url <- "https://npclassifier.gnps2.org"
order <- "/classify?smiles="
new <- new |>
  dplyr::mutate(query = RCurl::curlEscape(urls = structure_smiles_2D))
queries <- new$query
cached <- "&cached" # actually return wrong results?

log_debug("Classifying ...")
if (length(queries) != 0) {
  xs <- seq_len(length(queries))

  list_df <- getClass(xs = xs)

  df_new <- dplyr::bind_rows(list_df) |>
    dplyr::mutate_all(as.character)

  if (exists("old")) {
    df <- dplyr::bind_rows(old, df_new) |>
      dplyr::distinct()
  } else {
    df <- df_new
  }
} else {
  df <- old
}

taxonomy_semiclean <- treat_npclassifier_json()

df_semiclean_1 <- df |>
  dplyr::filter(!is.na(class)) |>
  dplyr::select(-pathway, -superclass) |>
  dplyr::full_join(taxonomy_semiclean, by = c("class" = "class"))

df_semiclean_2 <- df |>
  dplyr::filter(is.na(class) & !is.na(superclass)) |>
  dplyr::select(-pathway) |>
  dplyr::full_join(taxonomy_semiclean, by = c("superclass" = "superclass")) |>
  dplyr::select(-class.y) |>
  dplyr::rename(class = class.x)

df_semiclean_3 <- df |>
  dplyr::filter(is.na(class) & is.na(superclass))

df_semiclean <-
  rbind(df_semiclean_1, df_semiclean_2, df_semiclean_3) |>
  dplyr::filter(!is.na(query)) |>
  dplyr::filter(!is.na(class) |
    !is.na(superclass) |
    !is.na(pathway)) |>
  dplyr::distinct()

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimDictionariesStructureDictionaryNpclassifierFile)

log_debug("Exporting")
readr::write_delim(
  x = df_semiclean,
  delim = "\t",
  file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
