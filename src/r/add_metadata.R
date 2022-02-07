library(data.table)
library(DBI)
library(dplyr)
library(purrr)
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
    fread(file = pathDataInterimDictionariesStructureMetadata) %>%
    distinct(
      structure_inchikey = structureCleanedInchikey,
      structure_inchi = structureCleanedInchi,
      structure_smiles = structureCleanedSmiles,
      structure_molecular_formula = structureCleaned_molecularFormula,
      structure_exact_mass = structureCleaned_exactMass,
      structure_smiles_2D = structureCleaned_smiles2D,
      structure_nameIupac = structureCleaned_nameIupac,
      structure_nameTraditional = structureCleaned_nameTraditional,
      structure_stereocenters_total = structureCleaned_stereocenters_total,
      structure_stereocenters_unspecified = structureCleaned_stereocenters_unspecified
    ) %>%
    distinct(structure_inchikey,
      .keep_all = TRUE
    )

  log_debug(
    "We have",
    nrow(chemical_metadata),
    "metadata for structures"
  )

  chemical_taxonomy_1 <<-
    fread(file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile)

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
  ) %>%
    mutate_all(as.character)

  otl <- dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_otl"
  ) %>%
    mutate_all(as.character)

  meta <- dbGetQuery(
    conn = db,
    statement = "SELECT * FROM taxa_meta"
  ) %>%
    mutate_all(as.character)

  biological_metadata <- left_join(names, otl) %>%
    left_join(., meta, by = c("ott_id" = "id")) %>%
    filter(
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
    ) %>%
    distinct() %>%
    map_df(rev) %>%
    ## feeling it is better that way
    distinct(canonical_name, ott_id, rank, .keep_all = TRUE) %>%
    ## canonical_name important for synonyms
    pivot_wider(
      names_from = "rank",
      values_from = c("name", "unique_name.y", "ott_id.y")
    ) %>%
    select(
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
    ) %>%
    map_df(rev) %>%
    coalesce()

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
  df_complete <- df %>%
    left_join(., biological_metadata) %>%
    left_join(., chemical_metadata) %>%
    left_join(., chemical_taxonomy_1) %>%
    left_join(., chemical_taxonomy_2) %>%
    select(any_of(
      c(
        "structure_wikidata",
        "structure_inchikey",
        "structure_inchi",
        "structure_smiles",
        "structure_molecular_formula",
        "structure_exact_mass",
        "structure_smiles_2D",
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
