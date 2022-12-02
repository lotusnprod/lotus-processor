# title: "PAMDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readxl)
library(splitstackshape)
library(stringr)

# get paths
database <- databases$get("pamdb")

## files
data_original <-
  readxl::read_excel(path = database$sourceFiles$tsv) |>
  dplyr::mutate_all(as.character)

# selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid = MetID,
    structure_name = Name,
    structure_inchi = InChI,
    structure_smiles = SMILES,
    cas = `CAS number`,
    reference = References
  ) |>
  dplyr::mutate(organism_clean = "Pseudomonas aeruginosa")

data_manipulated <- data_selected |>
  splitstackshape::cSplit("reference",
    sep = "Pubmed:",
    fixed = TRUE,
    stripWhite = FALSE
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate(
    reference_title = str_extract(string = reference_1, pattern = "\".*\""),
    reference_original = ifelse(
      test = !is.na(reference_title),
      yes = NA,
      no = reference_1
    )
  ) |>
  splitstackshape::cSplit("reference_2",
    sep = " ",
    fixed = TRUE,
    stripWhite = FALSE
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::select(
    uniqueid,
    organism_clean,
    structure_name,
    structure_inchi,
    structure_smiles,
    cas,
    reference_original,
    reference_title,
    reference_pubmed = reference_2_02
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "pamdb",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_original",
      "reference_pubmed",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
