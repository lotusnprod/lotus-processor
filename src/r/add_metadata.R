library(DBI)
library(dplyr)
library(purrr)
library(readr)
library(RSQLite)

source("r/treat_npclassifier_taxonomy.R")
source("temp_classyfireTaxonomy.R")

#' Title
#'
#' @return
#' @export
#'
#' @examples
add_metadata <- function(df) {
  chemical_metadata <-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureMetadata,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    ) |>
    dplyr::distinct(
      structure_inchikey = structureCleanedInchikey,
      structure_inchi = structureCleanedInchi,
      structure_smiles = structureCleanedSmiles,
      structure_molecular_formula = structureCleaned_molecularFormula,
      structure_exact_mass = structureCleaned_exactMass,
      structure_xlogp = structureCleaned_xlogp,
      structure_smiles_2D = structureCleaned_smiles2D,
      structure_nameIupac = structureCleaned_nameIupac,
      structure_nameTraditional = structureCleaned_nameTraditional,
      structure_cid = structureCleaned_cid,
      structure_stereocenters_total = structureCleaned_stereocenters_total,
      structure_stereocenters_unspecified = structureCleaned_stereocenters_unspecified
    ) |>
    dplyr::distinct(structure_inchikey,
      structure_cid,
      .keep_all = TRUE
    )

  biological_metadata_2 <-
    readr::read_delim(
      file = pathDataInterimDictionariesOrganismMetadata,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    ) |>
    dplyr::filter(
      organismCleaned_dbTaxo == "GBIF Backbone Taxonomy" |
        organismCleaned_dbTaxo == "NCBI"
    ) |>
    dplyr::distinct(
      organism_name = organismCleaned,
      organismCleaned_id,
      organismCleaned_dbTaxo
    ) |>
    tidyr::pivot_wider(names_from = organismCleaned_dbTaxo, values_from = organismCleaned_id) |>
    dplyr::distinct(
      organism_name,
      organism_taxonomy_gbifid = `GBIF Backbone Taxonomy`,
      organism_taxonomy_ncbiid = NCBI
    ) |>
    dplyr::mutate_all(as.character)
  biological_metadata_2[biological_metadata_2 == "NULL"] <- NA

  log_debug(
    "We have",
    nrow(chemical_metadata),
    "metadata for structures"
  )

  chemical_taxonomy_1 <<-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )

  chemical_taxonomy_1 <<- treat_npclassifier_taxonomy()

  log_debug(
    "We have",
    nrow(chemical_taxonomy_1),
    "npclassifier classifications for structures"
  )

  log_debug(
    "We have",
    nrow(chemical_taxonomy_2),
    "classyfire classifications for structures"
  )

  drv <- SQLite()

  db <- dbConnect(
    drv = drv,
    dbname = pathDataInterimDictionariesOrganismDictionaryOTL
  )

  names <- dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_names"
  ) |>
    dplyr::mutate_all(as.character)

  otl <- dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_otl"
  ) |>
    dplyr::mutate_all(as.character)

  meta <- dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_meta"
  ) |>
    dplyr::mutate_all(as.character)

  biological_metadata <- dplyr::left_join(names, otl) |>
    dplyr::left_join(meta, by = c("ott_id" = "id")) |>
    dplyr::filter(
      rank %in% c(
        "domain",
        "kingdom",
        "phylum",
        "class",
        "order",
        "infraorder",
        "family",
        "subfamily",
        "tribe",
        "subtribe",
        "genus",
        "subgenus",
        "species",
        "subspecies",
        "varietas"
      )
    ) |>
    dplyr::distinct() |>
    purrr::map_df(rev) |>
    ## feeling it is better that way
    dplyr::distinct(canonical_name, ott_id, rank, .keep_all = TRUE) |>
    ## canonical_name important for synonyms
    tidyr::pivot_wider(
      names_from = "rank",
      values_from = c("name", "unique_name.y", "ott_id.y")
    ) |>
    dplyr::select(
      organism_name = canonical_name,
      organism_taxonomy_ottid = ott_id,
      organism_taxonomy_01domain = dplyr::matches("name_domain"),
      organism_taxonomy_02kingdom = dplyr::matches("name_kingdom"),
      organism_taxonomy_03phylum = dplyr::matches("name_phylum"),
      organism_taxonomy_04class = dplyr::matches("name_class"),
      organism_taxonomy_05order = dplyr::matches("name_order"),
      organism_taxonomy_06family = dplyr::matches("name_family"),
      organism_taxonomy_07tribe = dplyr::matches("name_tribe"),
      organism_taxonomy_08genus = dplyr::matches("name_genus"),
      organism_taxonomy_09species = dplyr::matches("name_species"),
      organism_taxonomy_10varietas = dplyr::matches("name_varietas")
    ) |>
    purrr::map_df(rev) |>
    dplyr::coalesce()

  if (nrow(biological_metadata) != 0) {
    biological_metadata[dplyr::setdiff(
      x = c(
        "organism_name",
        "organism_taxonomy_ottid",
        "organism_taxonomy_01domain",
        "organism_taxonomy_02kingdom",
        "organism_taxonomy_03phylum",
        "organism_taxonomy_04class",
        "organism_taxonomy_05order",
        "organism_taxonomy_06family",
        "organism_taxonomy_07tribe",
        "organism_taxonomy_08genus",
        "organism_taxonomy_09species",
        "organism_taxonomy_10varietas"
      ),
      y = names(biological_metadata)
    )] <- NA
  }

  log_debug(
    "We have",
    nrow(biological_metadata),
    "Open Tree of Life classifications for organisms"
  )

  log_debug("Adding useful metadata")
  df_complete <- df |>
    dplyr::left_join(biological_metadata) |>
    dplyr::left_join(biological_metadata_2) |>
    dplyr::left_join(chemical_metadata) |>
    dplyr::left_join(chemical_taxonomy_1) |>
    dplyr::left_join(chemical_taxonomy_2) |>
    dplyr::select(dplyr::any_of(
      c(
        "structure_wikidata",
        "structure_inchikey",
        "structure_inchi",
        "structure_smiles",
        "structure_molecular_formula",
        "structure_exact_mass",
        "structure_xlogp",
        "structure_smiles_2D",
        "structure_cid",
        "structure_nameIupac",
        "structure_nameTraditional",
        "structure_stereocenters_total",
        "structure_stereocenters_unspecified",
        "structure_taxonomy_npclassifier_01pathway",
        "structure_taxonomy_npclassifier_02superclass",
        "structure_taxonomy_npclassifier_03class",
        "structure_taxonomy_classyfire_chemontid",
        "structure_taxonomy_classyfire_01kingdom",
        "structure_taxonomy_classyfire_02superclass",
        "structure_taxonomy_classyfire_03class",
        "structure_taxonomy_classyfire_04directparent",
        "organism_wikidata",
        "organism_name",
        "organism_taxonomy_gbifid",
        "organism_taxonomy_ncbiid",
        "organism_taxonomy_ottid",
        "organism_taxonomy_01domain",
        "organism_taxonomy_02kingdom",
        "organism_taxonomy_03phylum",
        "organism_taxonomy_04class",
        "organism_taxonomy_05order",
        "organism_taxonomy_06family",
        "organism_taxonomy_07tribe",
        "organism_taxonomy_08genus",
        "organism_taxonomy_09species",
        "organism_taxonomy_10varietas",
        "reference_wikidata",
        "reference_doi",
        "manual_validation"
      )
    ))

  return(df_complete)
}
