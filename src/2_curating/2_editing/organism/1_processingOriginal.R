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

cut <-
  1000 # else verification takes too long and gnverifier times out

log_debug("  Step 1")
log_debug("... taxa ranks dictionary")
taxaRanksDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaRanks,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

wrongVerifiedDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaWrongVerified,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  as.list()

organismTable <-
  readr::read_delim(
    file = pathDataInterimTablesOriginalOrganismFull,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct()

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimTablesProcessedOrganism)
create_dir_with_rm(export = pathDataInterimTablesProcessedOrganismOriginal)

log_debug("submitting to GNVerifier")
if (.Platform$OS.type == "unix") {
  system(command = paste("bash", pathOriginalGnverifierScript))
} else {
  shell(paste("bash", pathOriginalGnverifierScript))
}

verified <-
  jsonlite::stream_in(con = file(
    pathDataInterimTablesProcessedOrganismVerifiedOriginalTable
  ))

verified_df <- verified |>
  data.frame() |>
  dplyr::select(-curation, -matchType) |>
  tidyr::unnest(results, names_repair = "minimal") |>
  dplyr::filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") |>
  dplyr::filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) |>
  dplyr::filter(isSynonym == FALSE) |>
  dplyr::mutate(organismType = "clean") |>
  dplyr::arrange(dplyr::desc(sortScore)) |>
  dplyr::distinct(name, dataSourceTitleShort, .keep_all = TRUE) |>
  dplyr::select(
    organismType,
    organismValue = name,
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

dataOrganismVerified <- dplyr::left_join(
  organismTable,
  verified_df
) |>
  dplyr::select(
    organismType,
    organismValue,
    organismCleaned,
    organismDbTaxo,
    taxonId,
    organismCleanedCurrent,
    taxonomy,
    rank
  )

dataOrganismNoVerified <- dataOrganismVerified |>
  dplyr::arrange(organismDbTaxo) |>
  dplyr::distinct(organismValue, .keep_all = TRUE) |>
  dplyr::filter(is.na(organismDbTaxo)) |>
  dplyr::distinct(organismValue) |>
  data.table()

dataOrganismVerified <- dataOrganismVerified |>
  dplyr::filter(!is.na(organismDbTaxo))

log_debug(pathDataInterimTablesOriginalOrganism)

if (nrow(dataOrganismNoVerified) != 0) {
  split_data_table(
    x = dataOrganismNoVerified,
    no_rows_per_frame = cut,
    # else verification takes too long and gnverifier times out
    text = "",
    path_to_store = pathDataInterimTablesOriginalOrganism
  )
}

log_debug("submitting to GNFinder")
if (.Platform$OS.type == "unix") {
  system(command = paste("bash", pathOriginalGnfinderScript))
} else {
  shell(paste("bash", pathOriginalGnfinderScript))
}

log_debug("treating GNFinder results")
length <-
  length(list.files(
    path = pathDataInterimTablesOriginalOrganism,
    pattern = "^[0-9]{6}.tsv"
  ))

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
    dplyr::bind_rows(dataCleanOriginalOrganism) |>
    dplyr::select(
      organismValue,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      dplyr::everything()
    ) |>
    dplyr::select(-ids, -dbQuality)
}

if (length(dataCleanOriginalOrganism) == 0) {
  dataCleanedOriginalOrganism <- data.frame() |>
    dplyr::mutate(
      organismValue = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxo = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
    )
}

dataCleanedOriginalOrganismUnique <- dataCleanedOriginalOrganism |>
  dplyr::distinct(organismValue, organismCleaned, .keep_all = TRUE)

log_debug("exporting ...")
if (length != 0) {
  log_debug(pathDataInterimTablesProcessedOrganismOriginalTable)
}

if (length != 0) {
  readr::write_delim(
    x = dataCleanedOriginalOrganism,
    delim = "\t",
    file = pathDataInterimTablesProcessedOrganismOriginalTable,
    na = ""
  )
}

if (length != 0) {
  log_debug(pathDataInterimTablesProcessedOrganismOriginalUniqueTable)
}

if (length != 0) {
  readr::write_delim(
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

readr::write_delim(
  x = dataOrganismVerified,
  delim = "\t",
  file = pathDataInterimTablesProcessedOrganismOriginalVerifiedTable,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
