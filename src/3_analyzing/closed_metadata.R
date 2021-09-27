source("r/log_debug.R")
log_debug("This script formats dnp as per ssot_union_wikidata.")
log_debug("It currently needs 'temp_classyfireTaxonomy.R' to be run before.")

start <- Sys.time()

safety <- FALSE

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")
source("r/y_as_na.R")
source("r/treat_npclassifier_taxonomy.R")
source("temp_classyfireTaxonomy.R")

library(data.table)
library(DBI)
library(dplyr)
library(purrr)
library(RSQLite)
library(tidyr)

log_debug("importing ...")
closed_pairs <-
  fread(file = pathDataInterimTablesAnalyzedClosedDbTriplets) %>%
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

wikidata_pairs <-
  fread(
    file = file.path(
      pathDataExternalDbSource,
      pathLastWdExport
    ),
    quote = ""
  ) %>%
  filter(!is.na(structure_inchikey) &
    !is.na(taxon_name) &
    !is.na(reference_doi)) %>%
  distinct(
    structure_wikidata = structure,
    structure_inchikey,
    organism_wikidata = taxon,
    organism_name = taxon_name,
    reference_wikidata = reference,
    reference_doi,
  )

log_debug(
  "We have",
  nrow(wikidata_pairs),
  "unique inchikey-taxon-doi vaidated triplets in wikidata"
)

chemical_metadata <-
  fread(file = pathDataInterimDictionariesStructureMetadata) %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    structure_inchi = structureCleanedInchi,
    structure_smiles = structureCleanedSmiles,
    structure_molecular_formula = structureCleaned_molecularFormula,
    structure_exact_mass = structureCleaned_exactMass,
    structure_smiles_2D = structureCleaned_smiles2D,
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

chemical_taxonomy_1 <-
  fread(file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile)

chemical_taxonomy_1 <- treat_npclassifier_taxonomy()

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
  pivot_wider(
    names_from = "rank",
    values_from = c("name", "unique_name.y", "ott_id.y")
  ) %>%
  select(
    organism_name = canonical_name,
    organism_taxonomy_ottid = ott_id,
    organism_taxonomy_01domain = name_domain,
    organism_taxonomy_02kingdom = name_kingdom,
    organism_taxonomy_03phylum = name_phylum,
    organism_taxonomy_04class = name_class,
    organism_taxonomy_05order = name_order,
    organism_taxonomy_06family = name_family,
    organism_taxonomy_07tribe = name_tribe,
    organism_taxonomy_08genus = name_genus,
    organism_taxonomy_09species = name_species,
    organism_taxonomy_10varietas = name_varietas
  ) %>%
  map_df(rev) %>%
  coalesce()

log_debug(
  "We have",
  nrow(biological_metadata),
  "Open Tree of Life classifications for organisms"
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
  "unique inchikey-taxon-doi vaidated triplets in dnp"
)

closed_only <-
  anti_join(closed_pairs, wikidata_pairs) %>% distinct()

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
  left_join(., biological_metadata) %>%
  left_join(., chemical_metadata) %>%
  left_join(., chemical_taxonomy_1) %>%
  left_join(., chemical_taxonomy_2) %>%
  select(
    structure_inchikey,
    structure_inchi,
    structure_smiles,
    structure_molecular_formula,
    structure_exact_mass,
    structure_smiles_2D,
    structure_stereocenters_total,
    structure_stereocenters_unspecified,
    structure_taxonomy_npclassifier_01pathway,
    structure_taxonomy_npclassifier_02superclass,
    structure_taxonomy_npclassifier_03class,
    structure_taxonomy_classyfire_chemontid,
    structure_taxonomy_classyfire_01kingdom,
    structure_taxonomy_classyfire_02superclass,
    structure_taxonomy_classyfire_03class,
    structure_taxonomy_classyfire_04directparent,
    organism_name,
    organism_taxonomy_ottid,
    organism_taxonomy_01domain,
    organism_taxonomy_02kingdom,
    organism_taxonomy_03phylum,
    organism_taxonomy_04class,
    organism_taxonomy_05order,
    organism_taxonomy_06family,
    organism_taxonomy_07tribe,
    organism_taxonomy_08genus,
    organism_taxonomy_09species,
    organism_taxonomy_10varietas
  )

if (safety == TRUE) {
  log_debug(
    "Exporting to",
    file.path(
      pathDataProcessed,
      pathLastFrozenClosed
    )
  )

  fwrite(
    x = closed_complete,
    file = file.path(
      pathDataProcessed,
      pathLastFrozenClosed
    )
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
