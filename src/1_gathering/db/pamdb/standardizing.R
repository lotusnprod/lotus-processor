# title: "PAMDB cleaneR"

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
database <- databases$get("pamdb")

##files
data_original <-
  read_excel(database$sourceFiles$tsv) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = MetID,
    name = Name,
    inchi = InChI,
    smiles = SMILES,
    cas = `CAS number`,
    reference = References
  ) %>%
  mutate(biologicalsource = "Pseudomonas aeruginosa")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "pam_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
database$writeInterim(data_standard)
