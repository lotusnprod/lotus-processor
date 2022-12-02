# title: "AntiBase cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)

# get paths
database <- databases$get("antibase")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv))

data_smiles <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$smi))

# selecting
data_selected <- data_original |>
  dplyr::select(
    MOL_ID,
    name = NAME,
    biologicalsource = SOURCE,
    reference_publishingDetails = REFERENCES
  ) |>
  dplyr::mutate(
    biologicalsource = gsub(
      pattern = "\\[+[[:alpha:]]+\\]",
      replacement = "",
      x = biologicalsource
    )
  ) |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = ";",
    direction = "long",
    fixed = TRUE
  ) |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = ",",
    direction = "long",
    fixed = TRUE
  )

data_corrected <- data_selected |>
  dplyr::left_join(data_smiles) |>
  dplyr::select(
    structure_name = name,
    organism_clean = biologicalsource,
    organism_dirty = biologicalsource,
    reference_publishingDetails,
    structure_smiles = SMILES
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected,
    db = "antibase",
    structure_field = "structure_smiles",
    organism_field = "organism_dirty",
    reference_field = "reference_publishingDetails"
  )

# exporting
database$writeInterim(data_standard)
