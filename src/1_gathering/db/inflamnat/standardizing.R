# title: "INFLAMNAT cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(readxl)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("inflamnat")

## files
data_original <-
  read_excel(path = database$sourceFiles$tsv,
             sheet = 1) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = Index,
    name = Name,
    smiles = SMILES,
    biologicalsource = Origin,
    pubchem = CID,
    reference = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "inf_1",
    structure_field = c("name", "smiles")
  )

# exporting
database$writeInterim(data_standard)

