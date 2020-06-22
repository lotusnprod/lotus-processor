# title: "treating bio"

# loading paths
source("paths.R")

## functions
source("functions/bio.R")
source("functions/helpers.R")
source("functions/parallel.R")
source("functions/log.R")

## libraries
library(data.table)
library(dplyr)
library(pbmcapply)
library(readr)
library(tidyverse)
library(tidyr)

log_debug("  Step 2")

## tcm names
tcmNamesDic <- read_delim(
  file = gzfile(pathDataInterimDictionariesTcmNames),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)
log_debug("      Loaded TCM names")
## common names
commonNamesDic <- read_delim(
  file = gzfile(pathDataInterimDictionariesCommonNames),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)
log_debug("      Loaded Common Names")
## black listed strings
blacklistDictionary <- read_delim(
  file = pathDataInterimDictionariesCommonBlackDic,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(n = str_count(string = blackName)) %>%
  arrange(desc(n)) %>%
  select(-n)

log_debug("      Loaded blacklist")
dataInterimOrganism <- dataCleanedOriginalOrganism %>%
  mutate(
    organismInterim = stri_replace_all_regex(
      str =  organismOriginal,
      pattern = paste("\\b", blacklistDictionary$blackName, "\\b", sep = ""),
      replacement = "",
      case_insensitive = TRUE,
      vectorize_all = FALSE
    )
  ) %>%
  rowwise() %>%
  mutate(organismInterim = stri_replace_all_regex(
    str =  organismInterim,
    pattern = paste("\\b", organismCleaned, "\\b", sep = ""),
    replacement = "",
  ))

log_debug("      Cleaned up")
dataInterimOrganism$organismInterim <-
  gsub(
    pattern = ".",
    replacement = "",
    x = dataInterimOrganism$organismInterim,
    fixed = TRUE
  )

dataInterimOrganism$organismInterim <-
  y_as_na(x = dataInterimOrganism$organismInterim, y = "")

dataInterimOrganism$organismInterim <-
  trimws(dataInterimOrganism$organismInterim)

dataInterimOrganismToFill <- dataInterimOrganism %>%
  filter(!is.na(organismInterim))

# manipulating organisms names
## removing some unnecessary characters
dataInterimOrganismToFill$organismInterim <- gsub(
  pattern = "_",
  replacement = " ",
  x = dataInterimOrganismToFill$organismInterim,
  fixed = TRUE
)

dataInterimOrganismToFill$organismInterim <- gsub(
  pattern = " ",
  replacement = " ",
  x = dataInterimOrganismToFill$organismInterim,
  fixed = TRUE
)

## replacing names
tcmNamesDic2 <- tcmNamesDic %>%
  filter(canonicalName != newCanonicalName)

log_debug("      Ready to replace")

replaceCommonNames <- function(value) {
  stri_replace_all_fixed(
    str_trim(value),
    c(
      commonNamesDic$vernacularName,
      tcmNamesDic$vernacularName,
      tcmNamesDic2$canonicalName
    ),
    c(
      commonNamesDic$canonicalName,
      tcmNamesDic$canonicalName,
      tcmNamesDic2$newCanonicalName
    ),
    case_insensitive = TRUE,
    vectorize_all = FALSE
  )
}

dataInterimOrganismToFill$organismInterim <- mclapply(
  FUN = replaceCommonNames,
  X = dataInterimOrganismToFill$organismInterim,
  mc.set.seed = TRUE,
  mc.silent = FALSE,
  mc.cores = numCores,
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

log_debug("       Finished replacing the common names")
dataInterimOrganismToFill$organismInterim <-
  as.character(dataInterimOrganismToFill$organismInterim)

dataInterimOrganismToFill$organismInterim <-
  y_as_na(x = dataInterimOrganismToFill$organismInterim, y = "")

dataInterimOrganismToFill$organismInterim <-
  y_as_na(x = dataInterimOrganismToFill$organismInterim, y = "NA")

dataInterimOrganismToFillGnfinder <- dataInterimOrganismToFill %>%
  data.table() %>%
  mutate_all(as.character) %>%
  filter(!is.na(organismInterim)) %>%
  select(organismInterim)
log_debug("     Exporting")

# exporting
split_data_table(
  x = dataInterimOrganismToFillGnfinder,
  no_rows_per_frame = cut,
  text = "translatedOrganismGnfinderUntil_",
  path_to_store = pathTranslatedOrganismDistinct
)
