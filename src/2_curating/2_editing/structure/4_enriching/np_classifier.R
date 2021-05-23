source("r/log_debug.R")
log_debug("This script adds np-classifier classes to structures dictionary")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(data.table)
library(jsonlite)
library(pbmcapply)
library(tidyverse)
library(RCurl)

log_debug("... functions")
source("r/getClass.R")
source("r/vroom_safe.R")
source("r/treat_npclassifier_json.R")

log_debug("loading smiles ...")
smiles <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureStereoCounted) %>%
  distinct(structure_smiles_2D = smilesSanitizedFlat) %>%
  tibble()
# smiles <-
#   vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata) %>%
#   distinct(structure_smiles_2D = structureCleaned_smiles2D) %>%
#   tibble()

log_debug("loading npClassifier taxonomy ...")
taxonomy <- fromJSON(txt = list.files(
  path = file.path(
    pathDataExternal,
    "taxonomySource/structure/npclassifier/"
  ),
  pattern = ".json",
  full.names = TRUE
))

old <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionaryNpclassifierFile) %>%
  distinct() %>%
  tibble()

new <- anti_join(smiles, old)
# new <- smiles

url <- "https://npclassifier.ucsd.edu"
order <- "/classify?smiles="
new <- new %>%
  mutate(query = curlEscape(structure_smiles_2D))
queries <- new$query
cached <- "&cached" # actually return wrong results?

if (length(queries) != 0) {
  X <- (seq_along(queries))

  list_df <- invisible(
    pbmclapply(
      FUN = getClass,
      X = X,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = parallel::detectCores() - 2,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

  df_new <- bind_rows(list_df) %>%
    mutate_all(as.character)

  df <- bind_rows(old, df_new) %>%
    distinct()
} else {
  df <- old
}

taxonomy_semiclean <- treat_npclassifier_json()

df_semiclean <- df %>%
  select(-pathway, -superclass) %>%
  full_join(taxonomy_semiclean) # %>%
# distinct(smiles, class, glycoside, query, superclass, pathway)

vroom_write_safe(
  x = df_semiclean,
  path = pathDataInterimDictionariesStructureDictionaryNpclassifierFile
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
