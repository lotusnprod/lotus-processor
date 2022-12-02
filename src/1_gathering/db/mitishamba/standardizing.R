# title: "MITISHAMBA cleaneR"

# loading paths
source("paths.R")
source("r/capitalize.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(readr)

# get paths
database <- databases$get("mitishamba")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv)) |>
  dplyr::mutate_all(as.character)

# selecting
data_selected <- data_original |>
  dplyr::select(
    smiles,
    biologicalsource = plant_species,
    name = common_name,
    reference_original = authors
  ) |>
  dplyr::mutate(reference_test = sub(
    pattern = "\\([0-9]{4}\\)",
    replacement = "§",
    x = reference_original
  )) |>
  splitstackshape::cSplit("reference_test", sep = "§") |>
  dplyr::mutate(reference_test_2 = sub(
    pattern = "^\\.",
    replacement = "",
    x = reference_test_2
  )) |>
  dplyr::mutate(reference_test_2 = sub(
    pattern = "^\\,",
    replacement = "",
    x = reference_test_2
  )) |>
  dplyr::mutate(reference_test_2 = trimws(x = reference_test_2)) |>
  dplyr::select(
    structure_name = name,
    biologicalsource,
    structure_smiles = smiles,
    reference_authors = reference_test_1,
    reference_original,
    reference_split = reference_test_2
  ) |>
  data.frame()

data_corrected <- data_selected |>
  splitstackshape::cSplit("biologicalsource", sep = ",", direction = "long") |>
  dplyr::filter(grepl(pattern = "[A-Z]", x = biologicalsource)) |>
  dplyr::mutate(organism_clean = capitalize(tolower(biologicalsource))) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected,
    db = "mitishamba",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_original",
      "reference_authors",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)
