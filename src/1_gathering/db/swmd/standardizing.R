# title: "SWMD cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("swmd")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

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
