# title: "dianatDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(readxl)
library(tidyverse)

# get paths
database <- databases$get("dianatdb")

## files
data_original <-
  read_excel(
    path = database$sourceFiles$tsv,
    sheet = 1
  ) %>%
  mutate_all(as.character)

# manipulating
data_manipulated <- data_original %>%
  select(
    structure_name = Name,
    structure_smiles = SMILE,
    organism_clean = Plant,
    reference_original = Citation
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "dianatdb",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)
