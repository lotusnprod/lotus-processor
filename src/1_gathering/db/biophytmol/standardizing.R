# title: "Biophytmol cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("biophytmol")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid,
    name,
    smiles,
    biologicalsource,
    reference
  ) %>%
  cSplit("biologicalsource", "     ") %>%
  cSplit("reference", "§") %>%
  cSplit("reference_1", ",", direction = "long") %>%
  cSplit("reference_2", ".", fixed = TRUE) %>%
  mutate(reference_authors = gsub("[0-9]\\) ", "", reference_2_01)) %>%
  select(
    uniqueid,
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource_1,
    reference_authors,
    reference_pubmed = reference_1,
    reference_title = reference_2_02,
    reference_journal = reference_2_03
  ) %>%
  mutate_all(as.character) %>%
  tibble()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_2",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c(
      "reference_authors",
      "reference_journal",
      "reference_pubmed",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
