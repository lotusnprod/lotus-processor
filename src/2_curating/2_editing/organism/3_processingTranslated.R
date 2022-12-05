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
  readr::read_delim(
    file = pathDataInterimTablesOriginalOrganismFull,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... translated organisms")
dataInterimOrganismToFill <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismTranslatedInterim,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... cleaned original organisms")
dataCleanedOriginalOrganism <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismOriginalTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... verified original organisms")
dataVerifiedOriginalOrganism <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismOriginalVerifiedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug(" ... taxa ranks dictionary")
taxaRanksDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaRanks,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("ensuring directories exist")
create_dir_with_rm(export = pathDataInterimTablesProcessedOrganismTranslated)

if (length(list.files(path = pathDataInterimTablesTranslatedOrganism, pattern = "tsv")) != 0) {
  log_debug("submitting to GNFinder")
  if (.Platform$OS.type == "unix") {
    system(command = paste("bash", pathTranslatedGnfinderScript))
  } else {
    shell(paste("bash", pathTranslatedGnfinderScript))
  }
}

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
    dplyr::bind_rows(dataCleanTranslatedOrganism) |>
    dplyr::select(
      organismInterim,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      dplyr::everything()
    ) |>
    dplyr::select(-ids, -dbQuality)
}

if (length(dataCleanTranslatedOrganism) == 0) {
  dataCleanedTranslatedOrganism <- data.frame() |>
    dplyr::mutate(
      organismInterim = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxon = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
    ) |>
    dplyr::mutate_all(as.character)
}

#' temporary fix as gnfinder does not allow the verification we want anymore
library(tidyr)
source("r/y_as_na.R")
wrongVerifiedDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaWrongVerified,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  as.list()
dataCleanedOrganismVerify <- dataCleanedTranslatedOrganism |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::distinct(organismCleaned)

readr::write_delim(
  x = dataCleanedOrganismVerify,
  file = gzfile(
    description = pathDataInterimTablesProcessedOrganismVerifyTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  na = "",
  delim = "\t",
  quote = "none",
  escape = "double"
)
## because gnverifier does not parse quotes

log_debug("submitting to GNVerifier")
if (.Platform$OS.type == "unix") {
  system(command = paste("bash", pathGnverifierScript))
} else {
  shell(paste("bash", pathGnverifierScript))
}

verified <-
  jsonlite::stream_in(con = file(pathDataInterimTablesProcessedOrganismVerifiedTable))

if (nrow(dataCleanedOrganismVerify != 0)) {
  verified_df <- verified |>
    data.frame() |>
    dplyr::select(-curation, -matchType) |>
    tidyr::unnest(results, names_repair = "minimal") |>
    dplyr::filter(dataSourceTitleShort != "IRMNG (old)" &
      dataSourceTitleShort != "IPNI") |>
    dplyr::filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) |>
    dplyr::arrange(dplyr::desc(sortScore)) |>
    dplyr::distinct(name, dataSourceTitleShort, .keep_all = TRUE) |>
    dplyr::select(
      organismCleaned = name,
      organismDbTaxo = dataSourceTitleShort,
      taxonId = currentRecordId,
      currentName,
      currentCanonicalFull,
      taxonomy = classificationPath,
      rank = classificationRanks
    )

  ## example ID 165 empty, maybe fill later on
  verified_df$organismDbTaxo <-
    y_as_na(verified_df$organismDbTaxo, "")
}

if (nrow(dataInterimOrganismToFill) != 0) {
  dataCleanedTranslatedOrganism2join <-
    dataInterimOrganismToFill |>
    dplyr::filter(!is.na(organismInterim)) |>
    dplyr::distinct(organismValue, organismInterim) |>
    dplyr::mutate_all(as.character)
}

if (nrow(dataInterimOrganismToFill) == 0) {
  dataCleanedTranslatedOrganism2join <- data.frame() |>
    dplyr::mutate(
      organismValue = NA,
      organismInterim = NA
    ) |>
    mutate_all(as.character)
}

if (length != 0) {
  dataCleanedTranslatedOrganismFull <-
    dplyr::left_join(
      dataCleanedTranslatedOrganism2join,
      dataCleanedTranslatedOrganism |>
        dplyr::distinct(organismInterim, organismCleaned)
    ) |>
    dplyr::left_join(
      dataCleanedTranslatedOrganism |>
        dplyr::distinct(
          organismCleaned,
          organismCleanedCurrent,
          organismDbTaxo,
          taxonId,
          taxonomy,
          rank,
        )
    ) |>
    dplyr::select(-organismInterim) |>
    dplyr::distinct(organismValue,
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
    dplyr::bind_rows(
      dataVerifiedOriginalOrganism |>
        dplyr::select(-organismType) |>
        dplyr::mutate_all(as.character),
      dataCleanedOriginalOrganism |>
        dplyr::mutate_all(as.character)
    )
}

dataCleanedOrganism <- dataCleanedOrganism |>
  dplyr::distinct(organismValue,
    organismCleaned,
    taxonId,
    .keep_all = TRUE
  ) |>
  dplyr::group_by(organismValue) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(organismCleaned) |
    !n > 1) |>
  dplyr::select(-n) |>
  dplyr::distinct(
    organismValue,
    organismCleaned
  )

dataCleanedOrganism <- organismTable_full |>
  dplyr::left_join(dataCleanedOrganism) |>
  dplyr::distinct()

log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedOrganismTranslatedTable)
readr::write_delim(
  x = dataCleanedOrganism,
  delim = "\t",
  file = pathDataInterimTablesProcessedOrganismTranslatedTable,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
