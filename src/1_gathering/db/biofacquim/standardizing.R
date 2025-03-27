# title: "BIOFACQUIM cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("biofacquim")

## files
data_original <-
  readr::read_delim(
    file = gzfile(description = database$sourceFiles$tsv),
    delim = "\t"
  ) |>
  dplyr::mutate_all(as.character)

## selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid = ID,
    structure_name = Name,
    structure_smiles = SMILES,
    organism_clean = Specie,
    reference_journal = Journal,
    reference_doi = DOI,
    reference_publishingDetails = Reference
  ) |>
  dplyr::mutate(
    organism_clean = gsub(
      pattern = "_",
      replacement = " ",
      x = organism_clean,
      fixed = TRUE
    )
  )

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "biofacquim",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_doi",
      "reference_journal"
    )
  )

# exporting
database$writeInterim(data_standard)
