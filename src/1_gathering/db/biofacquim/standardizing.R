# title: "BIOFACQUIM cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("biofacquim")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    name = Name,
    smiles = SMILES,
    biologicalsource = Specie,
    reference = Reference
  )

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_1",
    structure_field = c("name", "smiles")
  )

# exporting
database$writeInterim(data_standard)
