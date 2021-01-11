# title: "NANPDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("nanpdb")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
) %>%
  mutate_all(as.character) %>%
  mutate(id = row.names(.))

# selecting
data_selected <- data_original %>%
  mutate(name = NA) %>%
  select(
    uniqueid = id,
    name,
    smiles = SMILES.,
    biologicalsource = Source.Species.Information,
    pubchem = PubChem.,
    reference = Reference.information,
    reference_authors = Authors.information
  ) %>%
  cSplit("biologicalsource",
    " Known use:",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  select(uniqueid,
    name,
    smiles,
    biologicalsource = biologicalsource_1,
    pubchem,
    reference,
    reference_authors
  ) %>%
  mutate(biologicalsource = gsub(
    "Source: ",
    "",
    biologicalsource
  )) %>%
  mutate_all(as.character) %>%
  cSplit("reference",
    "Title: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  cSplit("reference_2",
    " PubMed: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  cSplit("reference_1",
    " Reference: ",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  cSplit("reference_authors",
    "Author(s): ",
    stripWhite = FALSE,
    fixed = TRUE
  ) %>%
  cSplit("reference_authors_2",
    "Curator(s): ",
    stripWhite = FALSE,
    fixed = TRUE
  ) %>%
  mutate(
    reference_title = reference_2_1,
    reference_authors = reference_authors_2_1,
    reference_publishingDetails = reference_1_2
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "nan_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_title",
      "reference_authors",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
