source("r/log_debug.R")
log_debug("This script performs canonical name recognition on the translated organism field.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... functions")
source("r/gnfinder_cleaning.R")

log_debug("loading ...")
log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("  Step 3")
log_debug("... files ...")
log_debug("full")
organismTable_full <-
  read_delim(file = pathDataInterimTablesOriginalOrganismFull)

log_debug("... translated organisms")
dataInterimOrganismToFill <-
  read_delim(file = pathDataInterimTablesCleanedOrganismTranslatedInterim)

log_debug("... cleaned original organisms")
dataCleanedOriginalOrganism <-
  read_delim(file = pathDataInterimTablesCleanedOrganismOriginalTable)

log_debug("... verified original organisms")
dataVerifiedOriginalOrganism <-
  read_delim(file = pathDataInterimTablesCleanedOrganismOriginalVerifiedTable)

log_debug(" ... taxa ranks dictionary")
taxaRanksDictionary <-
  read_delim(file = pathDataInterimDictionariesTaxaRanks)

log_debug("ensuring directories exist")
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

log_debug("submitting to GNFinder")
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

log_debug("cleaning GNFinder output")
if (length != 0) {
  for (i in num) {
    j <- i / cut
    log_debug(paste("step", j, "of", length))
    tryCatch(
      {
        dataCleanTranslatedOrganism[[j]] <-
          gnfinder_cleaning(
            num = i,
            organismCol = "organismInterim"
          )
      },
      error = function(e) {
        log_debug("ERROR :", conditionMessage(e))
      }
    )
  }
}

log_debug("selecting and reordering")
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
    select(-ids, -dbQuality)
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
    distinct(organismValue, organismInterim) %>%
    mutate_all(as.character)
}

if (nrow(dataInterimOrganismToFill) == 0) {
  dataCleanedTranslatedOrganism2join <- data.frame() %>%
    mutate(
      organismValue = NA,
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
    distinct(organismValue,
      organismCleaned,
      taxonId,
      .keep_all = TRUE
    )
}

if (length != 0) {
  dataCleanedOrganism <-
    bind_rows(
      dataVerifiedOriginalOrganism %>% select(-organismType),
      dataCleanedOriginalOrganism,
      dataCleanedTranslatedOrganismFull
    )
}

if (length == 0) {
  dataCleanedOrganism <-
    bind_rows(
      dataVerifiedOriginalOrganism %>% select(-organismType),
      dataCleanedOriginalOrganism
    )
}

dataCleanedOrganism <- dataCleanedOrganism %>%
  distinct(organismValue,
    organismCleaned,
    taxonId,
    .keep_all = TRUE
  ) %>%
  group_by(organismValue) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(organismCleaned) |
    !n > 1) %>%
  select(-n) %>%
  distinct(
    organismValue,
    organismCleaned
  )

dataCleanedOrganism <- dataCleanedOrganism %>%
  left_join(organismTable_full, .) %>%
  distinct()

log_debug("exporting ...")
log_debug(pathDataInterimTablesCleanedOrganismTranslatedTable)
write_delim(
  x = dataCleanedOrganism,
  file = pathDataInterimTablesCleanedOrganismTranslatedTable
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
