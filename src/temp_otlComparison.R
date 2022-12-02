source("r/log_debug.R")
log_debug("This script compares OTL (Open Tree of Life) IDs obtained via gnverify and via rotl API")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("loading ...")
log_debug("... libraries")
library(DBI)
library(dplyr)
library(readr)
library(rotl)
library(RSQLite)

canonical_name_colname <- "organismCleaned"

fullOrganisms <-
  readr::read_delim(file = pathDataInterimDictionariesOrganismDictionary)

drv <- RSQLite::SQLite()

db <- DBI::dbConnect(
  drv = drv,
  dbname = pathDataInterimDictionariesOrganismDictionaryOTL
)

otlVersion <- DBI::dbGetQuery(
  conn = db,
  statement = "SELECT
  taxa_names.canonical_name AS organismCleaned,
  taxa_otl.ott_id AS organismCleanedId
  FROM taxa_names
  LEFT JOIN taxa_otl
  ON taxa_names.search_string = taxa_otl.search_string"
) |>
  dplyr::filter(organismCleaned %in% fullOrganisms$organismCleaned) |>
  dplyr::filter(!is.na(organismCleanedId)) |>
  dplyr::distinct()

log_debug(nrow(otlVersion), "for rotl API version")

gnverifyVersion <- fullOrganisms |>
  dplyr::filter(organismCleaned_dbTaxo == "Open Tree of Life") |>
  dplyr::filter(!is.na(organismCleaned_id)) |>
  dplyr::distinct(organismCleaned, organismCleaned_id) |>
  dplyr::mutate(organismCleanedId = as.integer(organismCleaned_id))

log_debug(nrow(gnverifyVersion), "for gnverify version")

diff <- dplyr::anti_join(otlVersion, gnverifyVersion)

log_debug("so it seems that rotl version has \n", nrow(diff), "more results")

log_debug("but then if we have a closer look and only go for distinct IDs...")
log_debug("just an example to show what is meant:")
log_debug("applying filter(organismCleanedId == \"65272\") to both tables")

otlStrepto <- otlVersion |>
  filter(organismCleanedId == "65272")

log_debug("rotl table")
otlStrepto

gnverifyStrepto <- gnverifyVersion |>
  dplyr::filter(organismCleanedId == "65272")

log_debug("gnverify table")
gnverifyStrepto

log_debug("so when keeping only one name per ID...")

otlVersionDistinct <- otlVersion |>
  dplyr::left_join(gnverifyVersion |>
    dplyr::mutate(isInBoth = "Y")) |>
  ## to keep the same synonym for both
  dplyr::arrange(isInBoth) |>
  dplyr::distinct(organismCleanedId, .keep_all = TRUE)
log_debug(nrow(otlVersionDistinct), "distinct IDs for rotl API version")

gnverifyVersionDistinct <- gnverifyVersion |>
  dplyr::distinct(organismCleanedId, .keep_all = TRUE)

log_debug(
  nrow(gnverifyVersionDistinct),
  "distinct IDs for gnverify version"
)

diffDistinct_1 <-
  dplyr::anti_join(otlVersionDistinct, gnverifyVersionDistinct)

diffDistinct_2 <-
  dplyr::anti_join(gnverifyVersionDistinct, otlVersionDistinct)

log_debug(
  "so it seems that rotl version has \n",
  nrow(diffDistinct_1),
  "names that are not in gnverify"
)
log_debug(
  "and that gnverify version has \n",
  nrow(diffDistinct_2),
  "names that are not in rotl"
)

log_debug("Conclusion: almost the same but very interesting for synonyms filtering")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
