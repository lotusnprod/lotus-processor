# title: "UNPD cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("unpd")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceUnpdIntegrated),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
## atomizing references
data_selected <- data_original %>%
  mutate(
    reference = gsub("(\\(\\d+).\\s", "| ", ref),
    reference = sub("\\| ", "", reference)
  ) %>%
  select(
    biologicalsource = ln_reduced,
    reference,
    inchi = InChI,
    smiles = SMILES
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "unp_1",
    structure_field = c("inchi", "name", "smiles")
  )

# exporting
database$writeInterim(data_standard)
