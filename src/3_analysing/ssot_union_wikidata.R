cat(
  "This script verifies what we uploaded on Wikidata \n",
  "and compements it with some metadata."
)
cat("It currently needs 'temp_classyfireTaxonomy.R' to be run before. \n")

start <- Sys.time()

safety <- FALSE

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")
source("r/y_as_na.R")

library(data.table)
library(DBI)
library(RSQLite)
library(tidyverse)

cat("importing ... \n")
platinum_pairs <-
  fread(file = pathDataInterimTablesAnalysedPlatinum) %>%
  filter(
    !is.na(structureCleanedInchikey) &
      !is.na(organismCleaned) &
      !is.na(referenceCleanedDoi)
  ) %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned,
    reference_doi = referenceCleanedDoi
  ) %>%
  tibble()

cat(
  "We have",
  nrow(platinum_pairs),
  "unique inchikey-taxon-doi vaidated triplets in our SSOT \n"
)

manually_validated_pairs <-
  fread(file = "../data/validation/manuallyValidated.tsv.gz") %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    organism_name = organismCleaned,
    reference_doi = referenceCleanedDoi
  ) %>%
  mutate(manual_validation = "Y") %>%
  tibble()

platinum_pairs <-
  left_join(platinum_pairs, manually_validated_pairs)

cat(
  "Of which",
  nrow(platinum_pairs %>% filter(manual_validation == "Y")),
  "manually validated. \n"
)

wikidata_pairs <-
  fread(file = file.path(
    pathDataExternalDbSource,
    "210223_wikidata_query.tsv"
  )) %>%
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
    reference_title
  ) %>%
  tibble()

cat(
  "We have",
  nrow(wikidata_pairs),
  "unique inchikey-taxon-doi vaidated triplets in wikidata \n"
)

chemical_metadata <-
  fread(file = pathDataInterimDictionariesStructureMetadata) %>%
  distinct(
    structure_inchikey = structureCleanedInchikey,
    structure_inchi = structureCleanedInchi,
    structure_smiles = structureCleanedSmiles,
    structure_molecular_formula = structureCleaned_molecularFormula,
    structure_exact_mass = structureCleaned_exactMass
  ) %>%
  distinct(structure_inchikey,
    .keep_all = TRUE
  ) %>%
  ## needed
  tibble()

cat(
  "We have",
  nrow(chemical_metadata),
  "metadata for structures \n"
)

chemical_taxonomy_1 <-
  fread(file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile) %>%
  filter(!is.na(pathway)) %>%
  filter(
    pathway != "Fatty acids" |
      superclass != "Fatty acyls" |
      class != "Halogenated hydrocarbons" |
      grepl(pattern = "Cl|Br|I|F", x = smiles) ## because actually returned instead of NA/null
  ) %>%
  distinct(
    structure_smiles = smiles,
    structure_taxonomy_npclassifier_01pathway = pathway,
    structure_taxonomy_npclassifier_02superclass = superclass,
    structure_taxonomy_npclassifier_03class = class
  ) %>%
  group_by(structure_smiles) %>%
  summarize(across(
    c(
      "structure_taxonomy_npclassifier_01pathway",
      "structure_taxonomy_npclassifier_02superclass",
      "structure_taxonomy_npclassifier_03class"
    ),
    ~ gsub(
      pattern = "\\bNA\\b",
      replacement = "",
      x = paste(.x, collapse = "|")
    )
  )) %>%
  ungroup() %>%
  mutate(across(
    everything(),
    ~ y_as_na(x = .x, y = "")
  )) %>%
  distinct() %>%
  tibble()

cat(
  "We have",
  nrow(chemical_taxonomy_1),
  "npclassifier classifications for structures \n"
)

chemical_taxonomy_2 <- classy_temp %>%
  mutate(direct_parent = apply(.[, 3:13], 1, function(x) {
    tail(na.omit(x), 1)
  })) %>%
  distinct(
    structure_inchikey = inchikey,
    structure_taxonomy_classyfire_chemontid = chemont_id,
    structure_taxonomy_classyfire_01kingdom = `01kingdom`,
    structure_taxonomy_classyfire_02superclass = `02superclass`,
    structure_taxonomy_classyfire_03class = `03class`,
    structure_taxonomy_classyfire_04directparent = direct_parent,
  ) %>%
  tibble()

cat(
  "We have",
  nrow(chemical_taxonomy_2),
  "classyfire classifications for structures \n"
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
  distinct(ott_id, rank, .keep_all = TRUE) %>%
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
  coalesce() %>%
  tibble()

cat(
  "We have",
  nrow(biological_metadata),
  "Open Tree of Life classifications for organisms \n"
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
    reference_title,
    manual_validation
  )

cat(
  "We have",
  nrow(platinum_u_wd),
  "unique inchikey-taxon-doi vaidated triplets in both platinum and wikidata \n"
)

platinum_only <-
  anti_join(platinum_pairs, wikidata_pairs) %>% distinct()

# cat("We have",
#     nrow(platinum_only),
#     "unique inchikey-taxon-doi vaidated triplets present only in platinum \n")

# wd_only <-
#   anti_join(wikidata_pairs, platinum_pairs) %>% distinct()

# cat("We have",
#     nrow(wd_only),
#     "unique inchikey-taxon-doi vaidated triplets present only in wikidata \n")

cat("Adding useful metadata")
platinum_u_wd_complete <- platinum_u_wd %>%
  left_join(., biological_metadata) %>%
  left_join(., chemical_metadata) %>%
  left_join(., chemical_taxonomy_1) %>%
  left_join(., chemical_taxonomy_2) %>%
  select(
    structure_wikidata,
    structure_inchikey,
    structure_inchi,
    structure_smiles,
    structure_molecular_formula,
    structure_exact_mass,
    structure_taxonomy_npclassifier_01pathway,
    structure_taxonomy_npclassifier_02superclass,
    structure_taxonomy_npclassifier_03class,
    structure_taxonomy_classyfire_chemontid,
    structure_taxonomy_classyfire_01kingdom,
    structure_taxonomy_classyfire_02superclass,
    structure_taxonomy_classyfire_03class,
    structure_taxonomy_classyfire_04directparent,
    organism_wikidata,
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
    organism_taxonomy_10varietas,
    reference_wikidata,
    reference_doi,
    reference_title,
    manual_validation
  )

if (safety == TRUE) {
  cat(
    "Exporting to",
    file.path(
      pathDataProcessed,
      "210223_frozen.csv.gz"
    )
  )

  fwrite(
    x = platinum_u_wd,
    file = file.path(
      pathDataProcessed,
      "210223_frozen.csv.gz"
    )
  )

  cat(
    "Exporting to",
    file.path(
      pathDataProcessed,
      "210223_frozen_metadata.csv.gz"
    )
  )

  fwrite(
    x = platinum_u_wd_complete,
    file = file.path(
      pathDataProcessed,
      "210223_frozen_metadata.csv.gz"
    )
  )
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
