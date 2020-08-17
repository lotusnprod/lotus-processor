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
library(jsonlite)
library(tidyverse)

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
  no = file.remove(
    list.files(path = pathDataInterimTablesCleanedOrganismTranslated,
               full.names = TRUE)
  ) &
    dir.create(
      pathDataInterimTablesCleanedOrganismTranslated,
      showWarnings = FALSE
    )
)

system(command = paste("bash", pathTranslatedGnfinderScript))

length <-
  length(list.files(path = pathDataInterimTablesTranslatedOrganism,
                    pattern = 'tsv'))

if (length != 0)
  num <- as.integer(seq(
    from = 1 * cut,
    to = length * cut,
    by = cut
  ))

dataCleanTranslatedOrganism <- list()

# cleaning GNFinder output
if (length != 0)
  for (i in num) {
    j <- i / cut
    print(paste("step", j, "of", length))
    tryCatch({
      dataCleanTranslatedOrganism[[j]] <-
        gnfinder_cleaning(num = i,
                          organismCol = "organismInterim")
    }, error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    })
  }

# selecting and reordering
if (length(dataCleanTranslatedOrganism) != 0)
  dataCleanedTranslatedOrganism <-
  bind_rows(dataCleanTranslatedOrganism) %>%
  select(
    organismInterim,
    organismCleaned = canonicalname,
    organismCleanedCurrent = canonicalnameCurrent,
    organismDbTaxo = dbTaxo,
    everything()
  ) %>%
  select(-nchar, -sum)

if (length(dataCleanTranslatedOrganism) == 0)
  dataCleanedTranslatedOrganism <- data.frame() %>%
  mutate(
    organismInterim = NA,
    organismCleaned = NA,
    organismCleanedCurrent = NA,
    organismDbTaxon = NA,
    value_min = NA,
    value_max = NA,
    taxonId = NA,
    taxonomy = NA,
    rank = NA,
    ids = NA,
    dbQuality = NA
  ) %>%
  mutate_all(as.character)

if (nrow(dataInterimOrganismToFill) != 0)
  dataCleanedTranslatedOrganism2join <-
  dataInterimOrganismToFill %>%
  mutate(organismInterim = ifelse(
    test = organismInterim == word(organismOriginal, 3),
    yes = NA,
    no = organismInterim
  )) %>% # this is to avoid too big family groups because of some "Asteraceae" etc making caluclations too big
  filter(!is.na(organismInterim)) %>%
  distinct(organismOriginal, organismInterim) %>%
  mutate_all(as.character)

if (nrow(dataInterimOrganismToFill) == 0)
  dataCleanedTranslatedOrganism2join <- data.frame() %>%
  mutate(organismOriginal = NA,
         organismInterim = NA) %>%
  mutate_all(as.character)

if (length != 0)
  dataCleanedTranslatedOrganismFull <-
  left_join(dataCleanedTranslatedOrganism2join,
            dataCleanedTranslatedOrganism) %>%
  select(-organismInterim) %>%
  distinct(organismOriginal,
           organismCleaned,
           taxonId,
           .keep_all = TRUE)

if (length != 0)
  dataCleanedOrganism <-
  rbind(dataCleanedOriginalOrganism,
        dataCleanedTranslatedOrganismFull)

if (length != 0)
  dataCleanedOrganism <- dataCleanedOrganism %>%
  distinct(organismOriginal,
           organismCleaned,
           taxonId,
           .keep_all = TRUE) %>%
  group_by(organismOriginal) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(organismCleaned) |
           !n > 1) %>%
  select(-n)

print("manipulating taxonomic levels")

if (length != 0)
  dataCleanedOrganismManipulated <-
  manipulating_taxo(dfsel = dataCleanedOrganism,
                    dic = taxaRanksDictionary)

if (length == 0)
  dataCleanedOrganismManipulated <- data.frame() %>%
  mutate(
    organismOriginal = NA,
    organismCleaned = NA,
    organismDbTaxo = NA,
    organsimDbTaxoQuality = NA,
    organismTaxonIds = NA,
    organismTaxonRanks = NA,
    organismTaxonomy = NA,
    organism_1_kingdom = NA,
    organism_2_phylum = NA,
    organism_3_class = NA,
    organism_4_order = NA,
    organism_5_family = NA,
    organism_6_genus = NA,
    organism_7_species = NA,
    organism_8_variety = NA
  )

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
