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
library(pbmcapply)
library(readr)
library(rotl)
library(RSQLite)

canonical_name_colname <- "organismCleaned"

dataCuratedOrganismAuto <-
  read_delim(file = pathDataInterimTablesProcessedOrganismFinal)

if (works_locally_only == FALSE) {
  triplesPostWikidata <-
    read_delim(file = wikidataLotusExporterDataOutputTriplesPath)

  organismsPostWikidata <-
    read_delim(file = wikidataLotusExporterDataOutputTaxaPath)

  postWikidata <- left_join(
    triplesPostWikidata %>% distinct(taxon),
    organismsPostWikidata %>% distinct(wikidataId, names_pipe_separated),
    by = c("taxon" = "wikidataId")
  ) %>%
    select(organismCleaned = names_pipe_separated)

  dataCuratedOrganismAuto <- dataCuratedOrganismAuto %>%
    bind_rows(., postWikidata) %>%
    distinct()
}

new_matched_names <- dataCuratedOrganismAuto %>%
  drop_na(!!as.name(canonical_name_colname)) %>%
  distinct(!!as.name(canonical_name_colname)) %>%
  mutate(search_string = tolower(organismCleaned)) %>%
  distinct(
    !!as.name(canonical_name_colname),
    search_string
  ) %>%
  select(
    canonical_name := !!as.name(canonical_name_colname),
    search_string
  ) %>%
  data.frame()

drv <- SQLite()

db <- dbConnect(
  drv = drv,
  dbname = pathDataInterimDictionariesOrganismDictionaryOTL
)

previously_matched_names <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_names"
)

previously_matched_otl <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_otl"
)

previously_matched_meta <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_meta"
)

# dbListObjects(db)

# dbListFields(db, "taxa_names")
# dbListFields(db, "taxa_otl")
# dbListFields(db, "taxa_meta")

new_matched_names <-
  anti_join(new_matched_names, previously_matched_names)

taxa_names <- new_matched_names %>%
  distinct(search_string)

dbWriteTable(
  conn = db,
  name = "taxa_names",
  value = new_matched_names,
  row.names = FALSE,
  append = TRUE
)

new_matched_names <- taxa_names$search_string

if (is_empty(new_matched_names) == FALSE) {
  new_matched_otl_exact <- tnrs_match_names(
    names = new_matched_names,
    do_approximate_matching = FALSE,
    include_suppressed = FALSE
  )

  taxa_approx <- new_matched_otl_exact %>%
    filter(is.na(unique_name)) %>%
    drop_na(search_string) %>%
    distinct(search_string)

  new_matched_otl_exact <- new_matched_otl_exact %>%
    filter(!is.na(unique_name))

  ## very doubtful quality
  if (length(taxa_approx != 0)) {
    new_matched_otl_approx <- tnrs_match_names(
      names = taxa_approx$search_string,
      do_approximate_matching = TRUE,
      include_suppressed = FALSE
    )
  }

  ## not joining it for now since results of fuzzy seem really bad
  new_matched_otl <-
    bind_rows(new_matched_otl_exact) %>% ## new_matched_otl_approx
    data.frame() ## loosing some comments with df conversion


  dbWriteTable(
    conn = db,
    name = "taxa_otl",
    value = new_matched_otl,
    row.names = FALSE,
    append = TRUE
  )

  new_ott_id <- new_matched_otl %>%
    distinct(ott_id)

  X <- seq_along(1:round(nrow(new_ott_id) / 100))

  ott_list <- list()
  for (i in X) {
    ott_list[[i]] <-
      new_ott_id$ott_id[(i * 100 - 99):(i * 100)][!is.na(new_ott_id$ott_id[(i *
        100 - 99):(i * 100)])]
  }

  get_otl_lineage <- function(X) {
    tryCatch({
      taxon_info <- taxonomy_taxon_info(
        ott_ids = ott_list[[X]],
        include_lineage = TRUE,
        include_terminal_descendants = TRUE
      )

      taxon_lineage <- taxon_info %>%
        tax_lineage()

      list_df <- list()

      for (i in seq_along(1:length(taxon_lineage))) {
        list_df[[i]] <- bind_rows(
          data.frame(
            id = ott_list[[X]][i],
            rank = taxon_info[[i]]$rank,
            name = taxon_info[[i]]$name,
            unique_name = taxon_info[[i]]$unique_name,
            ott_id = as.character(taxon_info[[i]]$ott_id)
          ),
          data.frame(id = ott_list[[X]][i], taxon_lineage[[i]])
        )
      }

      df <- bind_rows(list_df)
      return(df)
    })
  }

  new_matched_meta_list <-
    pbmclapply(
      FUN = get_otl_lineage,
      X = X,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )

  new_matched_meta <- bind_rows(new_matched_meta_list)

  dbWriteTable(
    conn = db,
    name = "taxa_meta",
    value = new_matched_meta,
    row.names = FALSE,
    append = TRUE
  )
}

dbDisconnect(db)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
