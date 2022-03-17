source("r/log_debug.R")
log_debug(
  "This script verifies what we uploaded on Wikidata",
  "and complements it with some metadata."
)
log_debug("It currently needs 'temp_classyfireTaxonomy.R' to be run before.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")
source("r/add_metadata.R")
source("r/y_as_na.R")

library(DBI)
library(dplyr)
library(purrr)
library(readr)
library(RSQLite)
library(splitstackshape)
library(tidyr)

log_debug("importing ...")
platinum_pairs_raw <-
  read_delim(
    file = pathDataInterimTablesAnalyzedPlatinum,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  filter(
    !is.na(structureCleanedInchikey) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  )

platinum_pairs <- platinum_pairs_raw %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned,
    reference_doi = referenceCleanedDoi
  )

log_debug(
  "We have",
  nrow(platinum_pairs),
  "unique inchikey-taxon-doi vaidated triplets in our SSOT"
)

manually_validated_pairs <-
  read_delim(
    file = "../data/validation/manuallyValidated.tsv.gz",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned,
    reference_doi = referenceCleanedDoi
  ) %>%
  mutate(manual_validation = "Y")

platinum_pairs <-
  left_join(platinum_pairs, manually_validated_pairs)

log_debug(
  "Of which",
  nrow(platinum_pairs %>% filter(manual_validation == "Y")),
  "manually validated."
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

platinum_u_wd <-
  inner_join(platinum_pairs, wikidata_pairs) %>%
  distinct(
    structure_wikidata,
    structure_inchikey,
    organism_wikidata,
    organism_name,
    reference_wikidata,
    reference_doi,
    manual_validation
  )

platinum_no_wd <-
  anti_join(platinum_pairs, wikidata_pairs) %>%
  distinct(
    structureCleanedInchikey = structure_inchikey,
    organismCleaned = organism_name,
    referenceCleanedDoi = reference_doi
  ) %>%
  left_join(platinum_pairs_raw)

log_debug(
  "We have",
  nrow(platinum_u_wd),
  "unique inchikey-taxon-doi vaidated triplets in both platinum and wikidata"
)

platinum_only <-
  anti_join(platinum_pairs, wikidata_pairs) %>%
  distinct()

platinum_u_wd_complete <- add_metadata(df = platinum_u_wd)

if (safety == TRUE) {
  log_debug(
    "Exporting to",
    file.path(
      pathDataProcessed,
      gsub(
        pattern = "_metadata",
        replacement = "",
        x = pathLastFrozen
      )
    )
  )

  write_delim(
    x = platinum_u_wd,
    delim = "\t",
    file = file.path(
      pathDataProcessed,
      gsub(
        pattern = "_metadata",
        replacement = "",
        x = pathLastFrozen
      )
    )
  )

  write_delim(
    x = platinum_no_wd,
    delim = "\t",
    file = pathDataInterimTablesAnalyzedPlatinumNew
  )

  log_debug(
    "Exporting to",
    file.path(
      pathDataProcessed,
      pathLastFrozen
    )
  )

  write_delim(
    x = platinum_u_wd_complete,
    delim = "\t",
    file = file.path(
      pathDataProcessed,
      pathLastFrozen
    )
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
