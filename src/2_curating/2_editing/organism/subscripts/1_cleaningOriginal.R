# title: "cleaning original organism"

# loading
## paths
source("paths.R")

## functions
source("functions/bio.R")
source("functions/helpers.R")
source("functions/log.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R") # shouldnt these path be in the path.R ???
source("2_curating/2_editing/organism/functions/gnfinder_cleaning.R") # shouldnt these path be in the path.R ???

## libraries
library(data.table)
library(dplyr)
library(jsonlite)
library(readr)
library(stringr)
library(tidyverse)
library(tidyr)

log_debug("  Step 1")
# writing path
## dictionaries
### taxa levels
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesCleaned),
  dir.create(pathDataInterimTablesCleaned),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesCleanedOrganism),
  dir.create(pathDataInterimTablesCleanedOrganism),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesCleanedOrganismOriginal),
  dir.create(pathDataInterimTablesCleanedOrganismOriginal),
  FALSE
)

system(command = paste("bash", pathOriginalGnfinderScript))

length <-
  length(list.files(path = pathDataInterimTablesOriginalOrganism,
                    pattern = 'tsv'))

cut <- 10000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

dataCleanOriginalOrganism <- list()

# cleaning GNFinder output
for (i in num) {
  j <- i / cut
  tryCatch({
    dataCleanOriginalOrganism[[j]] <-
      gnfinder_cleaning(num = i,
                        organismCol = "organismOriginal")
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# selecting and reordering
dataCleanedOriginalOrganism <-
  bind_rows(dataCleanOriginalOrganism) %>%
  select(
    organismOriginal,
    organismCleaned = canonicalname,
    organismDbTaxo = dbTaxo,
    everything()
  ) %>%
  select(-nchar, -sum)

dataCleanedOriginalOrganismUnique <-
  dataCleanedOriginalOrganism %>%
  distinct(organismOriginal, organismCleaned, .keep_all = TRUE)

# exporting
write.table(
  x = dataCleanedOriginalOrganism,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismOriginalTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

write.table(
  x = dataCleanedOriginalOrganismUnique,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismOriginalUniqueTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
