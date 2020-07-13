# title: "cleaning translated organism"

# loading
## paths
source("paths.R")

## functions
source("functions/bio.R")
source("functions/helpers.R")
source("functions/log.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R")
source("2_curating/2_editing/organism/functions/gnfinder_cleaning.R")

## libraries
library(data.table)
library(dplyr)
library(jsonlite)
library(readr)
library(stringr)
library(tidyverse)
library(tidyr)

log_debug("  Step 3")

## interim
dataInterimOrganismToFill <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedInterim),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## cleaned original names
dataCleanedOriginalOrganism <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismOriginalTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

### taxa levels
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesCleanedOrganismTranslated),
  dir.create(pathDataInterimTablesCleanedOrganismTranslated),
  FALSE
)

system(command = paste("bash", pathTranslatedGnfinderScript))

length <-
  length(list.files(path = pathDataInterimTablesTranslatedOrganism,
                    pattern = 'tsv'))

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

dataCleanTranslatedOrganism <- list()

# cleaning GNFinder output
for (i in num) {
  j <- i / cut
  tryCatch({
    dataCleanTranslatedOrganism[[j]] <-
      gnfinder_cleaning(num = i,
                        organismCol = "organismInterim")
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# selecting and reordering
dataCleanedTranslatedOrganism <-
  bind_rows(dataCleanTranslatedOrganism) %>%
  select(
    organismInterim,
    organismCleaned = canonicalname,
    organismDbTaxo = dbTaxo,
    everything()
  ) %>%
  select(-nchar, -sum)

dataCleanedTranslatedOrganism2join <-
  dataInterimOrganismToFill %>%
  select(organismOriginal, organismInterim) %>%
  mutate_all(as.character)

dataCleanedTranslatedOrganismFull <-
  left_join(dataCleanedTranslatedOrganism2join,
            dataCleanedTranslatedOrganism) %>%
  select(-organismInterim)

dataCleanedOrganism <-
  rbind(dataCleanedOriginalOrganism,
        dataCleanedTranslatedOrganismFull)

dataCleanedOrganism <- dataCleanedOrganism %>%
  distinct(organismOriginal,
           organismCleaned,
           .keep_all = TRUE) %>%
  group_by(organismOriginal) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(organismCleaned) |
           !n > 1) %>%
  select(-n)

dataCleanedOrganismManipulated <-
  manipulating_taxo(dfsel = dataCleanedOrganism,
                    dic = taxaRanksDictionary) %>%
  select(-rank, -taxonomy)

# exporting
write.table(
  x = dataCleanedOrganismManipulated,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismTranslatedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
