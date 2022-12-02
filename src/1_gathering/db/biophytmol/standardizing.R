# title: "Biophytmol cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)

# get paths
database <- databases$get("biophytmol")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv)) |>
  dplyr::mutate_all(as.character)

# selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid,
    name,
    smiles,
    biologicalsource,
    reference
  ) |>
  splitstackshape::cSplit("biologicalsource", sep = "     ") |>
  splitstackshape::cSplit("reference", sep = "§") |>
  splitstackshape::cSplit("reference_1", sep = ",", direction = "long") |>
  splitstackshape::cSplit("reference_2", sep = ".", fixed = TRUE) |>
  dplyr::mutate(reference_authors = gsub(
    pattern = "[0-9]\\) ",
    replacement = "",
    x = reference_2_01
  )) |>
  dplyr::select(
    uniqueid,
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource_1,
    reference_authors,
    reference_pubmed = reference_1,
    reference_title = reference_2_02,
    reference_journal = reference_2_03
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::tibble()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "biophytmol",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_authors",
      "reference_journal",
      "reference_pubmed",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
