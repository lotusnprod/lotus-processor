cat("This script performs canonical name recognition on the translated organism field. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/manipulating_taxo.R")
source("r/gnfinder_cleaning.R")

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(tidyverse)

log_debug("  Step 3")
cat("... files ... \n")
cat("... translated organisms \n")
dataInterimOrganismToFill <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedInterim),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

cat("... cleaned original organisms \n")
dataCleanedOriginalOrganism <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismOriginalTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

cat(" ... taxa ranks dictionary \n")
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesCleanedOrganismTranslated),
  yes = dir.create(pathDataInterimTablesCleanedOrganismTranslated),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesCleanedOrganismTranslated,
      full.names = TRUE
    )
  ) &
    dir.create(
      pathDataInterimTablesCleanedOrganismTranslated,
      showWarnings = FALSE
    )
)

cat("submitting to GNFinder \n")
system(command = paste("bash", pathTranslatedGnfinderScript))

length <-
  length(list.files(
    path = pathDataInterimTablesTranslatedOrganism,
    pattern = "tsv"
  ))

cut <- 10000

if (length != 0) {
  num <- as.integer(seq(
    from = 1 * cut,
    to = length * cut,
    by = cut
  ))
}

dataCleanTranslatedOrganism <- list()

cat("cleaning GNFinder output \n")
if (length != 0) {
  for (i in num) {
    j <- i / cut
    cat(paste("step", j, "of", length))
    tryCatch(
      {
        dataCleanTranslatedOrganism[[j]] <-
          gnfinder_cleaning(
            num = i,
            organismCol = "organismInterim"
          )
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }
}

cat("selecting and reordering \n")
if (length(dataCleanTranslatedOrganism) != 0) {
  dataCleanedTranslatedOrganism <-
    bind_rows(dataCleanTranslatedOrganism) %>%
    select(
      organismInterim,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      everything()
    ) %>%
    select(-nchar, -sum) %>%
    distinct()
}

if (length(dataCleanTranslatedOrganism) == 0) {
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
}

if (nrow(dataInterimOrganismToFill) != 0) {
  dataCleanedTranslatedOrganism2join <-
    dataInterimOrganismToFill %>%
    filter(!is.na(organismInterim)) %>%
    distinct(organismOriginal, organismInterim) %>%
    mutate_all(as.character)
}

if (nrow(dataInterimOrganismToFill) == 0) {
  dataCleanedTranslatedOrganism2join <- data.frame() %>%
    mutate(
      organismOriginal = NA,
      organismInterim = NA
    ) %>%
    mutate_all(as.character)
}

if (length != 0) {
  dataCleanedTranslatedOrganismFull <-
    left_join(
      dataCleanedTranslatedOrganism2join,
      dataCleanedTranslatedOrganism %>% distinct(organismInterim, organismCleaned)
    ) %>%
    left_join(
      .,
      dataCleanedTranslatedOrganism %>% distinct(
        organismCleaned,
        organismCleanedCurrent,
        organismDbTaxo,
        taxonId,
        taxonomy,
        rank,
        ids,
        dbQuality
      )
    ) %>%
    select(-organismInterim) %>%
    distinct(organismOriginal,
      organismCleaned,
      taxonId,
      .keep_all = TRUE
    )
}

if (length != 0) {
  dataCleanedOrganism <-
    bind_rows(
      dataCleanedOriginalOrganism,
      dataCleanedTranslatedOrganismFull
    )
}

if (length != 0) {
  dataCleanedOrganism <- dataCleanedOrganism %>%
    distinct(organismOriginal,
      organismCleaned,
      taxonId,
      .keep_all = TRUE
    ) %>%
    group_by(organismOriginal) %>%
    add_count() %>%
    ungroup() %>%
    filter(!is.na(organismCleaned) |
      !n > 1) %>%
    select(-n)
}

cat("manipulating taxonomic levels \n")
if (length != 0 &
  nrow(dataCleanedOrganism) != 0) {
  dataCleanedOrganismManipulated <-
    manipulating_taxo(
      dfsel = dataCleanedOrganism,
      dic = taxaRanksDictionary
    )
}

if (length == 0 |
  nrow(dataCleanedOrganism) == 0) {
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
}

cat("exporting ... \n")
cat(pathDataInterimTablesCleanedOrganismTranslatedTable, "\n")
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
