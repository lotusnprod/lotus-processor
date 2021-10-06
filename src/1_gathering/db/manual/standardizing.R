# title: "manual input cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("manual")

args <- commandArgs(trailingOnly = TRUE)

## file
data_manual <- read_delim(file = args[1], delim = "\t")

## selecting
data_selected <- data_manual %>%
  select(
    structure_smiles = smiles,
    organism_clean = organism,
    reference_doi = doi
  )
# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "manual",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = "reference_doi"
  )

# exporting
ifelse(
  test = !dir.exists("../data/interim"),
  yes = dir.create("../data/interim"),
  no = paste("../data/interim", "exists")
)
ifelse(
  test = !dir.exists("../data/interim/manual"),
  yes = dir.create("../data/interim/manual"),
  no = paste("../data/interim/manual", "exists")
)
database$writeInterim(data_standard)
