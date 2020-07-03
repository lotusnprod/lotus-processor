# title: "datawarrior cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("datawarrior")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# manipulating
data_manipulated <- data_original %>%
  select(
    name,
    smiles = Smiles,
    inchi = InChI,
    biologicalsource = remarks,
    reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "dat_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
database$writeInterim(data_standard)
