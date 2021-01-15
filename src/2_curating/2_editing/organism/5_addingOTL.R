cat("This script retrieves data from Open Tree of Life (OTL) \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... libraries \n")
library(DBI)
library(rotl)
library(RSQLite)
library(tidyverse)

canonical_name_colname <- "organismCleaned"

dataCuratedOrganismAuto <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismFinal)

## temporary, will be written properly if validated
triplesPostWikidata <-
  vroom_read_safe(path = "../../wikidataLotusExporter/data/output/compound_reference_taxon.tsv")

organismsPostWikidata <-
  vroom_read_safe(path = "../../wikidataLotusExporter/data/output/taxa.tsv")

postWikidata <- left_join(
  triplesPostWikidata %>% distinct(taxon),
  organismsPostWikidata %>% distinct(wikidataId, names_pipe_separated),
  by = c("taxon" = "wikidataId")
) %>%
  select(organismCleaned = names_pipe_separated)

new_matched_names <- dataCuratedOrganismAuto %>%
  bind_rows(., postWikidata) %>%
  drop_na(!!as.name(canonical_name_colname)) %>%
  distinct(!!as.name(canonical_name_colname)) %>%
  mutate(search_string = tolower(
    sub(
      pattern = "(\\w+\\s+\\w+).*",
      replacement = "\\1",
      x = organismCleaned
    )
  )) %>%
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
new_matched_otl_approx <- tnrs_match_names(
  names = taxa_approx$search_string,
  do_approximate_matching = TRUE,
  include_suppressed = FALSE
)

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

new_matched_meta_list <-
  taxonomy_taxon_info(
    ott_ids = new_ott_id$ott_id,
    include_lineage = TRUE
  ) %>%
  tax_lineage()

new_matched_meta <- bind_rows(new_matched_meta_list,
  .id = "id"
)

dbWriteTable(
  conn = db,
  name = "taxa_meta",
  value = new_matched_meta,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)
