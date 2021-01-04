# title: "BIOFACQUIM cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("biofacquim")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    name = Name,
    smiles = SMILES,
    biologicalsource = Specie,
    reference_journal = Journal,
    reference_doi = DOI,
    reference_publishingDetails = Reference
  ) %>%
  mutate(biologicalsource = gsub(
    pattern = "_",
    replacement = " ",
    x = biologicalsource,
    fixed = TRUE
  ))

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_doi",
      "reference_journal",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)