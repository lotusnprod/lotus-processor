# title: "SWMD cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

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
    name = 2,
    pubchem = 3,
    chemspider = 4,
    biologicalsource = 8,
    geo = 9,
    extraction = 10,
    smiles = 14,
    inchi = 15,
    reference_original = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "swm_1",
    structure_field = c("name", "smiles", "inchi"),
    reference_field = c("reference_original")
  )

# exporting
database$writeInterim(data_standard)
