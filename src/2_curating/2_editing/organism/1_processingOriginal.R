source("r/log_debug.R")
log_debug("This script performs canonical name recognition on the original organism field.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("...paths")
source("paths.R")

log_debug("... functions")
source("r/gnfinder_cleaning.R")
source("r/split_data_table.R")
source("r/y_as_na.R")

log_debug("loading ...")
log_debug("... libraries")
library(data.table)
library(dplyr)
library(readr)
library(tidyr)

log_debug("  Step 1")
log_debug("... taxa ranks dictionary")
taxaRanksDictionary <-
  read_delim(file = pathDataInterimDictionariesTaxaRanks)

wrongVerifiedDictionary <-
  read_delim(
    file = pathDataInterimDictionariesTaxaWrongVerified,
    delim = "\t"
  ) %>%
  as.list()

organismTable <-
  read_delim(
    file = pathDataInterimTablesOriginalOrganismFull,
    delim = "\t"
  ) %>%
  distinct()

log_debug("ensuring directories exist")
ifelse(
  test = !dir.exists(pathDataInterimTablesProcessed),
  yes = dir.create(pathDataInterimTablesProcessed),
  no = paste(pathDataInterimTablesProcessed, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesProcessedOrganism),
  yes = dir.create(pathDataInterimTablesProcessedOrganism),
  no = paste(pathDataInterimTablesProcessedOrganism, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesProcessedOrganismOriginal),
  yes = dir.create(pathDataInterimTablesProcessedOrganismOriginal),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesProcessedOrganismOriginal,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesProcessedOrganismOriginal,
      showWarnings = FALSE
    )
)

log_debug("submitting to GNVerifier")
system(command = paste("bash", pathOriginalGnverifierScript))

verified <-
  stream_in(con = file(pathDataInterimTablesProcessedOrganismVerifiedOriginalTable))

verified_df <- verified %>%
  data.frame() %>%
  select(
    -curation,
    -matchType
  ) %>%
  unnest(preferredResults, names_repair = "minimal") %>%
  filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") %>%
  filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) %>%
  mutate(organismType = "clean") %>%
  select(
    organismType,
    organismValue = input,
    organismCleaned = currentCanonicalFull,
    organismDbTaxo = dataSourceTitleShort,
    taxonId = currentRecordId,
    organismCleanedCurrent = currentName,
    taxonomy = classificationPath,
    rank = classificationRanks
  )

## example ID 165 empty, maybe fill later on
verified_df$organismDbTaxo <-
  y_as_na(verified_df$organismDbTaxo, "")

dataOrganismVerified <- left_join(
  organismTable,
  verified_df
) %>%
  select(
    organismType,
    organismValue,
    organismCleaned,
    organismDbTaxo,
    taxonId,
    organismCleanedCurrent,
    taxonomy,
    rank
  )

dataOrganismNoVerified <- dataOrganismVerified %>%
  arrange(organismDbTaxo) %>%
  distinct(organismValue, .keep_all = TRUE) %>%
  filter(is.na(organismDbTaxo)) %>%
  distinct(organismValue) %>%
  data.table()

dataOrganismVerified <- dataOrganismVerified %>%
  filter(!is.na(organismDbTaxo))

log_debug(pathDataInterimTablesOriginalOrganism)

if (nrow(dataOrganismNoVerified) != 0) {
  split_data_table(
    x = dataOrganismNoVerified,
    no_rows_per_frame = 2000,
    # else verification takes too long and gnverifier times out
    text = "",
    path_to_store = pathDataInterimTablesOriginalOrganism
  )
}

log_debug("submitting to GNFinder")
system(command = paste("bash", pathOriginalGnfinderScript))

log_debug("treating GNFinder results")
length <-
  length(list.files(
    path = pathDataInterimTablesOriginalOrganism,
    pattern = "^[0-9]{6}.tsv"
  ))

cut <-
  2000 # else verification takes too long and gnverifier times out

if (length != 0) {
  num <- as.integer(seq(
    from = 1 * cut,
    to = length * cut,
    by = cut
  ))
}

dataCleanOriginalOrganism <- list()

log_debug("cleaning GNFinder output")
if (length != 0) {
  for (i in num) {
    j <- i / cut
    log_debug(paste("step", j, "of", length))
    tryCatch(
      {
        dataCleanOriginalOrganism[[j]] <-
          gnfinder_cleaning(
            num = i,
            organismCol = "organismValue"
          )
      },
      error = function(e) {
        log_debug("ERROR :", conditionMessage(e))
      }
    )
  }
}

log_debug("selecting and reordering")
if (length(dataCleanOriginalOrganism) != 0) {
  dataCleanedOriginalOrganism <-
    bind_rows(dataCleanOriginalOrganism) %>%
    select(
      organismValue,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      everything()
    ) %>%
    select(-ids, -dbQuality)
}

if (length(dataCleanOriginalOrganism) == 0) {
  dataCleanedOriginalOrganism <- data.frame() %>%
    mutate(
      organismValue = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxo = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
    )
}

dataCleanedOriginalOrganismUnique <- dataCleanedOriginalOrganism %>%
  distinct(organismValue, organismCleaned, .keep_all = TRUE)

log_debug("exporting ...")
if (length != 0) {
  log_debug(pathDataInterimTablesProcessedOrganismOriginalTable)
}

if (length != 0) {
  write_delim(
    x = dataCleanedOriginalOrganism,
    delim = "\t",
    file = pathDataInterimTablesProcessedOrganismOriginalTable
  )
}

if (length != 0) {
  log_debug(
    pathDataInterimTablesProcessedOrganismOriginalUniqueTable
  )
}

if (length != 0) {
  write_delim(
    x = dataCleanedOriginalOrganismUnique,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimTablesProcessedOrganismOriginalUniqueTable,
      compression = 9,
      encoding = "UTF-8"
    ),
    quote = "none",
    escape = "double"
  )
  ## because of univocity parser settings
}

write_delim(
  x = dataOrganismVerified,
  delim = "\t",
  file = pathDataInterimTablesProcessedOrganismOriginalVerifiedTable
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
