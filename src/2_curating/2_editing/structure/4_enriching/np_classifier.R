cat("This script adds np-classifier classes to structures dictionary \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(pbmcapply)
library(tidyverse)
library(RCurl)

cat("... functions \n")
source("r/getClass.R")
source("r/vroom_safe.R")

cat("loading smiles ... \n")
smiles <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureSmiles) %>%
  distinct() %>%
  tibble()
# smiles <-
#   vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary) %>%
#   distinct(smiles = structureCleanedSmiles) %>%
#   tibble()

if (works_locally_only == FALSE) {
  triplesPostWikidata <-
    vroom_read_safe(path = wikidataLotusExporterDataOutputTriplesPath)

  structuresPostWikidata <-
    vroom_read_safe(path = wikidataLotusExporterDataOutputStructuresPath)

  postWikidata <- left_join(
    triplesPostWikidata %>% distinct(compound),
    structuresPostWikidata %>% distinct(wikidataId, isomericSmiles, canonicalSmiles),
    by = c("compound" = "wikidataId")
  ) %>%
    mutate(smiles = ifelse(
      test = !is.na(isomericSmiles),
      yes = isomericSmiles,
      no = canonicalSmiles
    )) %>%
    drop_na(smiles) %>%
    distinct(smiles)

  smiles <- smiles %>%
    bind_rows(., postWikidata)
}

old <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionaryNpclassifierFile) %>%
  distinct() %>%
  tibble()

new <- anti_join(smiles, old)
# new <- smiles

url <- "https://npclassifier.ucsd.edu"
order <- "/classify?smiles="
queries <- curlEscape(new$smiles)
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
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

  df_new <- bind_rows(list_df) %>% mutate_all(as.character)

  df <- bind_rows(old, df_new) %>%
    distinct()
} else {
  df <- old
}

vroom_write_safe(
  x = df,
  path = pathDataInterimDictionariesStructureDictionaryNpclassifierFile
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
