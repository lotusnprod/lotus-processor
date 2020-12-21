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

cat("loading smiles ... \n")
smiles <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedStructureSmiles),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  distinct() %>%
  tibble()

old <- read_delim(
  file = gzfile(
    pathDataInterimDictionariesStructureDictionaryNpclassifierFile
  ),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
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

  df_new <- bind_rows(list_df)

  df <- bind_rows(old, df_new) %>%
    distinct()
} else {
  df <- old
}

fwrite(
  df,
  file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile,
  quote = FALSE,
  sep = "\t",
  showProgress = TRUE
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
