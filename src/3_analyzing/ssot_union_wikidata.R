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
  readr::read_delim(
    file = pathDataInterimTablesAnalyzedPlatinum,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::filter(
    !is.na(structureCleanedInchikey) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  )

platinum_pairs <- platinum_pairs_raw |>
  dplyr::distinct(
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
  readr::read_delim(
    file = "../data/validation/manuallyValidated.tsv.gz",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned,
    reference_doi = referenceCleanedDoi
  ) |>
  dplyr::mutate(manual_validation = "Y")

platinum_pairs <-
  dplyr::left_join(platinum_pairs, manually_validated_pairs)

log_debug(
  "Of which",
  nrow(
    platinum_pairs |>
      dplyr::filter(manual_validation == "Y")
  ),
  "manually validated."
)

data_organism <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputTaxaPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "organism_wikidata" = "wikidataId",
      "organismCleaned" = "names_pipe_separated"
    )
  ) |>
  splitstackshape::cSplit("organismCleaned", sep = "|", direction = "long") |>
  dplyr::distinct()

data_structures <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputStructuresPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structure_wikidata" = "wikidataId",
      "structureCleanedInchiKey" = "inchiKey",
      "structureCleanedInchi" = "inchi",
      "canonicalSmiles",
      "isomericSmiles"
    )
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    structureValue = if_else(
      condition = !is.na(isomericSmiles),
      true = isomericSmiles,
      false = canonicalSmiles
    )
  )

data_structures_translation <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureDictionary
  )

# data_structures_1 <- data_structures |>
#   dplyr::filter(!grepl(pattern = "\\|", x = structureCleanedInchiKey))
#
# data_structures_2 <- data_structures |>
#   dplyr::filter(grepl(pattern = "\\|", x = structureCleanedInchiKey)) |>
#   splitstackshape::cSplit(c("structureCleanedInchiKey", "structureCleanedInchi"),
#     sep = "|"
#   ) |>
#   tidyr::pivot_longer(cols = contains(c(
#     "structureCleanedInchiKey", "structureCleanedInchi"
#   ))) |>
#   dplyr::filter(!is.na(value)) |>
#   splitstackshape::cSplit("name",
#     sep = "_"
#   ) |>
#   tidyr::pivot_wider(names_from = "name_1") |>
#   dplyr::select(-name_2) |>
#   dplyr::distinct()

# data_structures <- data_structures_1 |>
#   dplyr::bind_rows(data_structures_2)

data_references <-
  readr::read_delim(
    file = wikidataLotusExporterDataOutputReferencesPath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "reference_wikidata" = "wikidataId",
      "referenceCleanedDoi" = "dois_pipe_separated",
      "referenceCleanedTitle" = "title",
    )
  ) |>
  splitstackshape::cSplit(
    "referenceCleanedDoi",
    sep = "|",
    direction = "long"
  ) |>
  dplyr::distinct()

wikidata_pairs <-
  readr::read_delim(
    file = pathLastWdExport,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    col_select = c(
      "structureValue" = "structure_smiles",
      "organismCleaned" = "organism_clean",
      "referenceCleanedDoi" = "reference_doi"
    )
  ) |>
  filter(
    !is.na(structureValue) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  ) |>
  dplyr::left_join(data_structures_translation) |>
  dplyr::left_join(data_organism) |>
  dplyr::left_join(data_structures) |>
  dplyr::left_join(data_references) |>
  dplyr::distinct(
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
  dplyr::inner_join(platinum_pairs, wikidata_pairs) |>
  dplyr::distinct(
    structure_wikidata,
    structure_inchikey,
    organism_wikidata,
    organism_name,
    reference_wikidata,
    reference_doi,
    manual_validation
  )

platinum_no_wd <-
  dplyr::anti_join(platinum_pairs, wikidata_pairs) |>
  dplyr::distinct(
    structureCleanedInchikey = structure_inchikey,
    organismCleaned = organism_name,
    referenceCleanedDoi = reference_doi
  ) |>
  dplyr::left_join(platinum_pairs_raw) |>
  ## subspecies not pushed to WD yet
  dplyr::filter(
    !grepl(pattern = "subspecies", x = organismCleaned_dbTaxoTaxonRanks)
  ) |>
  dplyr::filter(organismCleaned %in% data_organism$organismCleaned)

log_debug(
  "We have",
  nrow(platinum_u_wd),
  "unique inchikey-taxon-doi vaidated triplets in both platinum and wikidata"
)

platinum_only <-
  dplyr::anti_join(platinum_pairs, wikidata_pairs) |>
  dplyr::distinct()

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

  readr::write_delim(
    x = platinum_u_wd,
    delim = ",",
    file = file.path(
      pathDataProcessed,
      gsub(
        pattern = "_metadata",
        replacement = "",
        x = pathLastFrozen
      )
    )
  )

  readr::write_delim(
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

  readr::write_delim(
    x = platinum_u_wd_complete,
    delim = ",",
    file = file.path(
      pathDataProcessed,
      pathLastFrozen
    )
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
