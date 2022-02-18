source("r/log_debug.R")
log_debug("This script computes trees for pubchem from the frozen dataset.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("...paths")
source("paths.R")

library(data.tree)
library(dplyr)
library(networkD3)
library(readr)
library(tidyr)

LIMIT <- 999999
options(max.print = LIMIT)

lotus <- readr::read_delim(file = file.path(
  pathDataProcessed,
  pathLastFrozen
))

#' since PubChem matching is based on NCBI taxonomy
organisms_ids <- readr::read_delim(file = pathDataInterimDictionariesOrganismMetadata) |>
  dplyr::filter(organismCleaned_dbTaxo == "NCBI") |>
  dplyr::distinct(organism_name = organismCleaned, organism_ncbi_id = organismCleaned_id)

classified <-
  readr::read_delim(file = pathDataInterimDictionariesStructureDictionaryNpclassifierFile)

structures_classified <- lotus |>
  dplyr::select(
    structure_id = structure_inchikey,
    # chemical_pathway = structure_taxonomy_npclassifier_01pathway,
    # chemical_superclass = structure_taxonomy_npclassifier_02superclass,
    # chemical_class = structure_taxonomy_npclassifier_03class,
    structure_smiles_2D
  ) |>
  dplyr::distinct()

organisms_classified <- lotus |>
  dplyr::select(
    organism_name,
    organism_01_domain = organism_taxonomy_01domain,
    organism_02_kingdom = organism_taxonomy_02kingdom,
    organism_03_phylum = organism_taxonomy_03phylum,
    organism_04_class = organism_taxonomy_04class,
    organism_05_order = organism_taxonomy_05order,
    organism_06_family = organism_taxonomy_06family,
    organism_07_tribe = organism_taxonomy_07tribe,
    organism_08_genus = organism_taxonomy_08genus,
    organism_09_species = organism_taxonomy_09species
  ) |>
  dplyr::distinct() |>
  dplyr::inner_join(organisms_ids)

#' on pubchem's team request: collapse not classified nodes
structures_cleaned <- structures_classified |>
  dplyr::left_join(classified) |>
  dplyr::distinct(
    structure_id,
    chemical_pathway = pathway,
    chemical_superclass = superclass,
    chemical_class = class
  ) |>
  dplyr::mutate(chemical_pathway = ifelse(
    test = is.na(chemical_pathway),
    yes = "Not classified",
    no = chemical_pathway
  )) |>
  dplyr::mutate(chemical_superclass = ifelse(
    test = is.na(chemical_superclass),
    yes = ifelse(
      test = chemical_pathway == "Not classified",
      yes = chemical_superclass,
      no = "Not classified"
    ),
    no = chemical_superclass
  )) |>
  dplyr::mutate(chemical_class = ifelse(
    test = is.na(chemical_class),
    yes = ifelse(
      test = chemical_superclass == "Not classified" |
        is.na(chemical_superclass),
      yes = chemical_class,
      no = "Not classified"
    ),
    no = chemical_class
  )) |>
  data.frame()

organisms_cleaned <- organisms_classified |>
  data.frame()

#' on pubchem's team request: remove not classified nodes
#' name of columns does not match anymore but not important
organisms_cleaned[] <-
  t(apply(organisms_cleaned, 1, function(x) {
    `length<-`(na.omit(x), length(x))
  }))

# organisms_cleaned[is.na(organisms_cleaned)] <- "Not classified"

structures_cleaned <- structures_cleaned |>
  dplyr::mutate(root = "Root") |>
  tidyr::unite(
    col = pathString,
    root,
    chemical_pathway,
    chemical_superclass,
    chemical_class,
    structure_id,
    sep = "§",
    na.rm = TRUE
  )

organisms_cleaned <- organisms_cleaned |>
  dplyr::mutate(root = "Root") |>
  tidyr::unite(
    col = pathString,
    root,
    organism_01_domain,
    organism_02_kingdom,
    organism_03_phylum,
    organism_04_class,
    organism_05_order,
    organism_06_family,
    organism_07_tribe,
    organism_08_genus,
    organism_09_species,
    organism_ncbi_id,
    sep = "§",
    na.rm = TRUE
  )

tree_chem_txt <-
  data.tree::as.Node(x = structures_cleaned, pathDelimiter = "§")

tree_chem_list <- data.tree::ToListExplicit(x = tree_chem_txt)[2]

tree_chem_json <- jsonlite::toJSON(tree_chem_list)

tree_bio_txt <-
  data.tree::as.Node(x = organisms_cleaned, pathDelimiter = "§")

tree_bio_list <- data.tree::ToListExplicit(x = tree_bio_txt)[2]

tree_bio_json <- jsonlite::toJSON(tree_bio_list)

sink(file = file.path(
  pathDataProcessed,
  pathLastTreeChemo
))
tree_chem_json
sink()

sink(file = file.path(
  pathDataProcessed,
  pathLastTreeBio
))
tree_bio_json
sink()
