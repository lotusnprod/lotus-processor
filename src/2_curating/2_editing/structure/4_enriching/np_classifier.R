cat("This script adds np-classifier classes to structures dictionary \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(pbmcapply)
library(tidyverse)

cat("... functions \n")
source("r/getClass.R")
source("r/vroom_safe.R")

cat("loading smiles ... \n")
smiles <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureSmiles) %>%
  distinct() %>%
  tibble()

old <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionaryNpclassifierFile) %>%
  distinct() %>%
  tibble()

new <- anti_join(smiles, old)

url <- "https://npclassifier.ucsd.edu"
order <- "/classify?smiles="
queries <- new$smiles
cached <- "&cached" # actually return wrong results?

if (length(queries) != 0) {
  X <- (1:length(queries))

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