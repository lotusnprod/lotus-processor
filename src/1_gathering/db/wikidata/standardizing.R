# title: "wikidata cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(tidyverse)
library(vroom)

# get paths
database <- databases$get("wikidata")

## files
data_organism <-
  vroom(
    file = wikidataLotusExporterDataOutputTaxaPath,
    delim = "\t"
  )

data_structures <-
  vroom(
    file = wikidataLotusExporterDataOutputStructuresPath,
    delim = "\t"
  )

data_references <-
  vroom(
    file = wikidataLotusExporterDataOutputReferencesPath,
    delim = "\t"
  )

data_triples <-
  vroom(
    file = wikidataLotusExporterDataOutputTriplesPath,
    delim = "\t"
  )

data_temp <-
  vroom(
    file = file.path(
      pathDataExternalDbSource,
      "210426_wikidata_query.tsv"
    ),
    delim = "\t"
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

data_manipulated <- data_temp %>%
  distinct(
    structure_inchi = compound_inchi,
    organism_clean = taxon_name,
    reference_doi
  )

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
