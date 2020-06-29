# title: "treating bio"

# loading
## paths
source("paths.R")

## functions
source("functions/bio.R")
source("functions/helpers.R")
source("functions/log.R")
source("2_curating/2_editing/bio/functions/manipulating_taxo.R") # shouldnt these path be in the path.R ???
source("2_curating/2_editing/bio/functions/gnfinder_cleaning.R") # shouldnt these path be in the path.R ???

## libraries
library(data.table)
library(dplyr)
library(jsonlite)
library(readr)
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


system(command = paste("bash", pathOriginalGnfinderScript))

length <- length(list.files(path = pathOriginalOrganismDistinct,
                            pattern = 'tsv'))

cut <- 10000

num <- as.integer(seq(®
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
