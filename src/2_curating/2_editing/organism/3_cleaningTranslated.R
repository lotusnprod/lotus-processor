cat("This script performs canonical name recognition on the translated organism field. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/gnfinder_cleaning.R")
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(tidyverse)

log_debug("  Step 3")
cat("... files ... \n")
cat("... translated organisms \n")
dataInterimOrganismToFill <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismTranslatedInterim)

cat("... cleaned original organisms \n")
dataCleanedOriginalOrganism <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismOriginalTable)

cat("... verified original organisms \n")
dataVerifiedOriginalOrganism <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismOriginalVerifiedTable)

cat(" ... taxa ranks dictionary \n")
taxaRanksDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesTaxaRanks)

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
    cat(paste("step", j, "of", length, "\n"))
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
    select(-nchar, -sum, -value_min, -value_max, -ids, -dbQuality)
}

if (length(dataCleanTranslatedOrganism) == 0) {
  dataCleanedTranslatedOrganism <- data.frame() %>%
    mutate(
      organismInterim = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxon = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
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
      dataVerifiedOriginalOrganism,
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
    select(-n) %>%
    distinct(
      organismOriginal,
      organismCleaned
    )
}

cat("exporting ... \n")
cat(pathDataInterimTablesCleanedOrganismTranslatedTable, "\n")
vroom_write_safe(
  x = dataCleanedOrganism,
  path = pathDataInterimTablesCleanedOrganismTranslatedTable
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
