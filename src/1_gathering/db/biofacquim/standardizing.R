# title: "BIOFACQUIM cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("biofacquim")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    structure_name = Name,
    structure_smiles = SMILES,
    organism_clean = Specie,
    reference_journal = Journal,
    reference_doi = DOI,
    reference_publishingDetails = Reference
  ) %>%
  mutate(organism_clean = gsub(
    pattern = "_",
    replacement = " ",
    x = organism_clean,
    fixed = TRUE
  ))

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "biofacquim",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c(
      "reference_doi",
      "reference_journal",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
