source("r/log_debug.R")
log_debug("This script helps finding nice examples")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("loading ...")
log_debug("... validated db, if running fullmode, this may take a while")
openDb <- readr::read_delim(
  file = gzfile(pathDataInterimTablesAnalyzedPlatinum),
  col_types = cols(.default = "c"),
  locale = locales,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

structureSearch_1 <- openDb |>
  dplyr::filter(structureType == "nominal") |>
  dplyr::distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch_2 <- openDb |>
  dplyr::filter(structureType == "smiles") |>
  dplyr::distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch_3 <- openDb |>
  dplyr::filter(structureType == "inchi") |>
  dplyr::distinct(structureValue,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

structureSearch <-
  rbind(structureSearch_1, structureSearch_2, structureSearch_3) |>
  dplyr::group_by(structureCleanedInchikey) |>
  dplyr::add_count() |>
  dplyr::arrange(dplyr::desc(n))

saltSearch <- structureSearch_3 |>
  dplyr::distinct(structureValue,
    .keep_all = TRUE
  ) |>
  dplyr::group_by(structureCleanedInchikey) |>
  dplyr::add_count() |>
  dplyr::filter(grepl(pattern = "\\.", x = structureValue)) |>
  dplyr::arrange(dplyr::desc(n))

maybeHit_salt <- openDb |>
  dplyr::filter(!is.na(referenceCleanedDoi)) |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::filter(structureCleanedInchikey == "KRKNYBCHXYNGOX-UHFFFAOYSA-N") |>
  dplyr::distinct(structureValue)

maybeHit_str <- openDb |>
  dplyr::filter(!is.na(referenceCleanedDoi)) |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N")

hitNames_str <- maybeHit_str |>
  dplyr::filter(structureType == "nominal") |>
  dplyr::distinct(structureValue)

hitSmiles_str <- maybeHit_str |>
  dplyr::filter(structureType == "smiles") |>
  dplyr::distinct(structureValue)

hitInchi_str <- maybeHit_str |>
  dplyr::filter(structureType == "inchi") |>
  dplyr::distinct(structureValue)

organismSearch <- openDb |>
  dplyr::distinct(organismType,
    organismValue,
    organismCleaned,
    .keep_all = TRUE
  ) |>
  dplyr::group_by(organismCleaned) |>
  dplyr::add_count() |>
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::filter(
    !grepl(pattern = "Streptomyces", x = organismCleaned) &
      !grepl(pattern = "Aspergillus", x = organismCleaned) &
      !grepl(pattern = "Fusarium", x = organismCleaned) &
      !grepl(pattern = ".*ae", x = organismCleaned)
  )

maybeHit_org <- openDb |>
  dplyr::filter(organismCleaned == "Oryza sativa")

hitNames_org <- maybeHit_org |>
  dplyr::distinct(
    organismType,
    organismValue
  )

reference <- openDb |>
  dplyr::filter(!is.na(referenceCleanedDoi))

referenceSearch <- reference |>
  dplyr::distinct(referenceValue,
    .keep_all = TRUE
  ) |>
  dplyr::group_by(referenceCleanedDoi) |>
  dplyr::add_count() |>
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::filter(
    !grepl(pattern = "Khimiya", x = referenceCleanedTitle) &
      !grepl(pattern = "Flavone and flavonol glycosides", x = referenceCleanedTitle) &
      !grepl(pattern = "The Handbook of Natural Flavonoids", x = referenceCleanedTitle) &
      n <= 10 &
      n >= 3
  )

maybeHit_ref <- openDb |>
  dplyr::filter(referenceCleanedDoi == "10.1021/np0600595")

hitNames_org <- maybeHit_org |>
  dplyr::distinct(
    organismType,
    organismValue
  )

doubleTest <- openDb |>
  dplyr::filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") |>
  dplyr::distinct(
    organismType,
    organismValue, organismCleaned
  ) |>
  dplyr::group_by(organismCleaned) |>
  dplyr::count() |>
  dplyr::arrange(dplyr::desc(n))

pairTest <- openDb |>
  dplyr::filter(structureCleanedInchikey == "OVSQVDMCBVZWGM-DTGCRPNFSA-N") |>
  dplyr::filter(organismCleaned == "Crataegus monogyna") |>
  dplyr::distinct(
    organismType,
    organismValue, structureValue, organismCleaned
  )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
