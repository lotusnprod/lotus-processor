source("r/log_debug.R")
log_debug("This script formats dnp as per ssot_union_wikidata.")
log_debug("It currently needs 'temp_classyfireTaxonomy.R' to be run before.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")
source("r/add_metadata.R")
source("r/y_as_na.R")
source("r/treat_npclassifier_taxonomy.R")
source("temp_classyfireTaxonomy.R")

library(DBI)
library(dplyr)
library(purrr)
library(readr)
library(RSQLite)
library(splitstackshape)
library(tidyr)

log_debug("importing ...")
closed_pairs <-
  read_delim(
    file = pathDataInterimTablesAnalyzedClosedDbTriplets,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  filter(!is.na(structureCleanedInchikey) &
    !is.na(organismCleaned)) %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned
  )

log_debug(
  "Closed DBs have",
  nrow(closed_pairs),
  "unique inchikey-taxon pairs"
)

data_organism <-
  read_delim(
    file = wikidataLotusExporterDataOutputTaxaPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "organism_wikidata" = "wikidataId",
      "organismCleaned" = "names_pipe_separated"
    )
  ) %>%
  cSplit("organismCleaned",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

data_structures <-
  read_delim(
    file = wikidataLotusExporterDataOutputStructuresPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structure_wikidata" = "wikidataId",
      "structureCleanedInchiKey" = "inchiKey",
      "structureCleanedInchi" = "inchi"
    )
  ) %>%
  distinct()

data_references <-
  read_delim(
    file = wikidataLotusExporterDataOutputReferencesPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "reference_wikidata" = "wikidataId",
      "referenceCleanedDoi" = "dois_pipe_separated",
      "referenceCleanedTitle" = "title",
    )
  ) %>%
  cSplit("referenceCleanedDoi",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

wikidata_pairs <-
  read_delim(
    file = pathLastWdExport,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structureCleanedInchi" = "structure_inchi",
      "organismCleaned" = "organism_clean",
      "referenceCleanedDoi" = "reference_doi"
    )
  ) %>%
  filter(
    !is.na(structureCleanedInchi) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  ) %>%
  left_join(., data_organism) %>%
  left_join(., data_structures) %>%
  left_join(., data_references) %>%
  distinct(
    structure_wikidata,
    structure_inchikey = structureCleanedInchiKey,
    organism_wikidata,
    organism_name = organismCleaned,
    reference_wikidata,
    reference_doi = referenceCleanedDoi,
  )

log_debug(
  "We have",
  nrow(wikidata_pairs),
  "unique inchikey-taxon-doi vaidated triplets in wikidata"
)

closed_u_wd <-
  left_join(closed_pairs, wikidata_pairs) %>%
  distinct(
    structure_wikidata,
    structure_inchikey,
    organism_wikidata,
    organism_name,
    reference_wikidata,
    reference_doi
  )

log_debug(
  "We have",
  nrow(closed_u_wd),
  "unique inchikey-taxon-doi vaidated triplets in closed resources"
)

closed_only <- anti_join(closed_pairs, wikidata_pairs) %>%
  distinct()

# log_debug("We have",
#     nrow(platinum_only),
#     "unique inchikey-taxon-doi vaidated triplets present only in platinum")

# wd_only <-
#   anti_join(wikidata_pairs, platinum_pairs) %>% distinct()

# log_debug("We have",
#     nrow(wd_only),
#     "unique inchikey-taxon-doi vaidated triplets present only in wikidata")

log_debug("Adding useful metadata")
closed_complete <- closed_only %>%
  add_metadata()

if (safety == TRUE) {
  log_debug(
    "Exporting to",
    file.path(
      pathDataProcessed,
      pathLastFrozenClosed
    )
  )

  write_delim(
    x = closed_complete,
    delim = ",",
    file = file.path(
      pathDataProcessed,
      pathLastFrozenClosed
    )
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
