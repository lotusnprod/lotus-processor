# title: "ANPDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(readr)

# get paths
database <- databases$get("anpdb")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv)) %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate(id = row.names(.))

# selecting
data_selected <- data_original |>
  dplyr::mutate(name = NA) |>
  dplyr::select(
    uniqueid = id,
    name,
    smiles = SMILES.,
    biologicalsource = Source.Species.Information,
    pubchem = PubChem.,
    reference = Reference.information,
    reference_authors = Authors.information
  ) |>
  splitstackshape::cSplit("biologicalsource",
    " Known use:",
    stripWhite = FALSE,
    fixed = FALSE
  ) |>
  dplyr::select(
    uniqueid,
    structure_name = name,
    structure_smiles = smiles,
    biologicalsource = biologicalsource_1,
    pubchem,
    reference,
    reference_authors
  ) |>
  dplyr::mutate(organism_clean = gsub(
    "Source: ",
    "",
    biologicalsource
  )) |>
  dplyr::mutate_all(as.character) |>
  splitstackshape::cSplit("reference",
    "Title: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) |>
  splitstackshape::cSplit("reference_2",
    " PubMed: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) |>
  splitstackshape::cSplit("reference_1",
    " Reference: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) |>
  splitstackshape::cSplit("reference_authors",
    "Author(s): ",
    stripWhite = FALSE,
    fixed = TRUE
  ) |>
  splitstackshape::cSplit("reference_authors_2",
    "Curator(s): ",
    stripWhite = FALSE,
    fixed = TRUE
  ) |>
  dplyr::mutate(
    reference_title = reference_2_1,
    reference_authors = reference_authors_2_1,
    reference_publishingDetails = reference_1_2
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "nanpdb",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_title",
      "reference_authors",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
