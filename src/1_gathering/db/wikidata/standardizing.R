# title: "wikidata cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)

# get paths
database <- databases$get("wikidata")

## files
data_organism <-
  read_delim(
    file = wikidataLotusExporterDataOutputTaxaPath
  ) %>%
  cSplit("names_pipe_separated",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

data_structures <-
  read_delim(
    file = wikidataLotusExporterDataOutputStructuresPath
  )

data_references <-
  read_delim(
    file = wikidataLotusExporterDataOutputReferencesPath
  ) %>%
  cSplit("dois_pipe_separated",
    sep = "|",
    direction = "long"
  ) %>%
  distinct()

data_triples <-
  read_delim(
    file = wikidataLotusExporterDataOutputTriplesPath
  )

# manipulating
data_manipulated <- data_triples %>%
  left_join(data_organism, by = c("taxon" = "wikidataId")) %>%
  left_join(data_structures, by = c("compound" = "wikidataId")) %>%
  left_join(data_references, by = c("reference" = "wikidataId")) %>%
  mutate(structure_smiles = if_else(
    condition = !is.na(isomericSmiles),
    true = isomericSmiles,
    false = canonicalSmiles
  )) %>%
  select(
    structure_smiles,
    structure_inchi = inchi,
    organism_clean = names_pipe_separated,
    reference_doi = dois_pipe_separated,
    reference_title = title
  ) %>%
  distinct() %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "wikidata",
    structure_field = "structure_inchi",
    organism_field = "organism_clean",
    reference_field = "reference_doi"
  )

# exporting
database$writeInterim(data_standard)
