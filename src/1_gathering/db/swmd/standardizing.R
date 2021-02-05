# title: "SWMD cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("swmd")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = 1,
    structure_name = 2,
    pubchem = 3,
    chemspider = 4,
    organism_clean = 8,
    geo = 9,
    extraction = 10,
    structure_smiles = 14,
    structure_inchi = 15,
    reference_original = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "swmd",
    structure_field = c("structure_name", "structure_smiles", "structure_inchi"),
    organism_field = "organism_clean",
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)
