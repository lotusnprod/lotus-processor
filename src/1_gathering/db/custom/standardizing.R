# title: "custom input cleaneR"

mode_custom <- TRUE

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("custom")

args <- commandArgs(trailingOnly = TRUE)

## file
data_custom <- readr::read_delim(file = args[1], delim = "\t")

## selecting
data_selected <- data_custom |>
  dplyr::select(
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = organism,
    reference_doi = doi
  )
# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "custom",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = "reference_doi"
  )

# exporting
create_dir(export = file.path(pathDataInterim, "custom"))
database$writeInterim(data_standard)
