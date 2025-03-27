# title: "STREPTOMEDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)

# get paths
database <- databases$get("streptomedb")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv))

## selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid,
    name,
    smiles,
    pubchem,
    reference_pubmed = pubmedid,
    biologicalsource
  ) |>
  splitstackshape::cSplit("reference_pubmed", sep = ";", direction = "long") |>
  splitstackshape::cSplit("biologicalsource", sep = ";", direction = "long") |>
  data.frame() |>
  dplyr::select(
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "streptomedb",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = "reference_pubmed"
  )

# exporting
database$writeInterim(data_standard)
