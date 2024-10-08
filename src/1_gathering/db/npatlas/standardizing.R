# title: "NPATLAS cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("npatlas")

## files
data_original <- readr::read_delim(file = database$sourceFiles$tsv)

# selecting
data_selected <- data_original |>
  dplyr::mutate(biologicalsource = paste(genus, origin_species, sep = " ")) |>
  dplyr::select(
    npaid,
    structure_name = compound_names,
    structure_inchi = compound_inchi,
    structure_smiles = compound_smiles,
    organism_clean = biologicalsource,
    reference_authors = original_reference_author_list,
    reference_doi = original_reference_doi,
    reference_journal = original_journal_title,
    reference_pubmed = original_reference_pmid,
    reference_title = original_reference_title
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npatlas",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_authors",
      "reference_doi",
      "reference_journal",
      "reference_pubmed",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
