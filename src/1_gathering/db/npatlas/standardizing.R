# title: "NPATLAS cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("npatlas")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  mutate(
    organism = paste(Genus,
                     `Origin Species`,
                     sep = " "),
    reference_doi = `Isolation Reference DOI`,
    reference_original = `Isolation Reference Citation`
  ) %>%
  select(
    NPAID,
    name = Names,
    InChIKey,
    inchi = InChI,
    smiles = SMILES,
    biologicalsource = organism,
    reference_doi,
    reference_original
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_2",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c("reference_doi", "reference_original")
  )

# exporting
database$writeInterim(data_standard)
