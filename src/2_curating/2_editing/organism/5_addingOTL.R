source("r/log_debug.R")
log_debug("This script retrieves data from Open Tree of Life (OTL)")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("loading ...")
log_debug("... libraries")
library(DBI)
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(purrr)
library(readr)
library(rotl)
library(RSQLite)
library(tidyr)

source("r/progressr.R")

canonical_name_colname <- "organismCleaned"

dataCuratedOrganismAuto <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismFinal,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

if (works_locally_only == FALSE) {
  triplesPostWikidata <-
    readr::read_delim(
      file = wikidataLotusExporterDataOutputTriplesPath,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )

  organismsPostWikidata <-
    readr::read_delim(
      file = wikidataLotusExporterDataOutputTaxaPath,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )

  postWikidata <- dplyr::left_join(
    triplesPostWikidata |>
      dplyr::distinct(taxon),
    organismsPostWikidata |>
      dplyr::distinct(wikidataId, names_pipe_separated),
    by = c("taxon" = "wikidataId")
  ) |>
    dplyr::select(organismCleaned = names_pipe_separated)

  dataCuratedOrganismAuto <- dataCuratedOrganismAuto |>
    dplyr::bind_rows(postWikidata) |>
    dplyr::distinct()
}

new_matched_names <- dataCuratedOrganismAuto |>
  tidyr::drop_na(!!as.name(canonical_name_colname)) |>
  dplyr::distinct(!!as.name(canonical_name_colname)) |>
  dplyr::mutate(search_string = tolower(organismCleaned)) |>
  dplyr::distinct(
    !!as.name(canonical_name_colname),
    search_string
  ) |>
  dplyr::select(
    canonical_name := !!as.name(canonical_name_colname),
    search_string
  ) |>
  data.frame()

drv <- RSQLite::SQLite()

create_dir(export = pathDataInterimDictionariesOrganismDictionaryOTL)

db <- DBI::dbConnect(
  drv = drv,
  dbname = pathDataInterimDictionariesOrganismDictionaryOTL
)

if ("taxa_names" %in% dbListTables(db)) {
  previously_matched_names <- DBI::dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_names"
  )

  new_matched_names <-
    dplyr::anti_join(new_matched_names, previously_matched_names)
}

if ("taxa_otl" %in% dbListTables(db)) {
  previously_matched_otl <- DBI::dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_otl"
  )
}

if ("taxa_meta" %in% dbListTables(db)) {
  previously_matched_meta <- DBI::dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_meta"
  )
}

taxa_names <- new_matched_names |>
  dplyr::distinct(search_string)

if ("taxa_names" %in% dbListTables(db)) {
  DBI::dbAppendTable(
    conn = db,
    name = "taxa_names",
    value = new_matched_names,
  )
} else {
  DBI::dbCreateTable(
    conn = db,
    name = "taxa_names",
    fields = new_matched_names,
  )
}

new_matched_names <- taxa_names$search_string

if (is_empty(new_matched_names) == FALSE) {
  new_matched_otl_exact <- rotl::tnrs_match_names(
    names = new_matched_names,
    do_approximate_matching = FALSE,
    include_suppressed = FALSE
  )

  taxa_approx <- new_matched_otl_exact |>
    dplyr::filter(is.na(unique_name)) |>
    tidyr::drop_na(search_string) |>
    dplyr::distinct(search_string)

  new_matched_otl_exact <- new_matched_otl_exact |>
    dplyr::filter(!is.na(unique_name))

  ## very doubtful quality
  if (length(taxa_approx != 0)) {
    new_matched_otl_approx <- rotl::tnrs_match_names(
      names = taxa_approx$search_string,
      do_approximate_matching = TRUE,
      include_suppressed = FALSE
    )
  }

  ## not joining it for now since results of fuzzy seem really bad
  new_matched_otl <-
    dplyr::bind_rows(new_matched_otl_exact) |> ## new_matched_otl_approx
    data.frame() ## loosing some comments with df conversion

  if ("taxa_otl" %in% dbListTables(db)) {
    DBI::dbAppendTable(
      conn = db,
      name = "taxa_otl",
      value = new_matched_otl
    )
  } else {
    DBI::dbCreateTable(
      conn = db,
      name = "taxa_otl",
      fields = new_matched_otl
    )
  }

  new_ott_id <- new_matched_otl |>
    dplyr::distinct(ott_id)

  xs <- seq_along(1:round(nrow(new_ott_id) / 100))

  ott_list <- list()
  for (i in xs) {
    ott_list[[i]] <-
      new_ott_id$ott_id[(i * 100 - 99):(i * 100)][!is.na(new_ott_id$ott_id[(i *
        100 - 99):(i * 100)])]
  }

  get_otl_lineage <- function(xs) {
    p <- progressr::progressor(along = xs)
    future.apply::future_lapply(
      future.seed = TRUE,
      X = xs,
      FUN = function(x) {
        tryCatch({
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          taxon_info <- rotl::taxonomy_taxon_info(
            ott_ids = ott_list[[x]],
            include_lineage = TRUE,
            include_terminal_descendants = TRUE
          )

          taxon_lineage <- taxon_info |>
            rotl::tax_lineage()

          list_df <- list()

          for (i in seq_along(1:length(taxon_lineage))) {
            list_df[[i]] <- dplyr::bind_rows(
              data.frame(
                id = ott_list[[x]][i],
                rank = taxon_info[[i]]$rank,
                name = taxon_info[[i]]$name,
                unique_name = taxon_info[[i]]$unique_name,
                ott_id = as.character(taxon_info[[i]]$ott_id)
              ),
              data.frame(id = ott_list[[x]][i], taxon_lineage[[i]])
            )
          }

          df <- dplyr::bind_rows(list_df)
          return(df)
        })
      }
    )
  }

  new_matched_meta_list <-
    get_otl_lineage(xs = xs) |>
    progressr::with_progress()

  new_matched_meta <- dplyr::bind_rows(new_matched_meta_list)

  DBI::dbWriteTable(
    conn = db,
    name = "taxa_meta",
    value = new_matched_meta,
    row.names = FALSE,
    append = TRUE
  )
}

DBI::dbDisconnect(db)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
