# title: "TRIFORC cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("triforc")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv2,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# cleaning
data_selected <- data_original %>%
  select(
    name = Name,
    cas = CAS,
    pubchem = `PubChem CID`,
    biologicalsource = Plant
  ) %>%
  cSplit("biologicalsource", ", ") %>%
  pivot_longer(
    4:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "biologicalsource",
    values_drop_na = TRUE
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tri_1",
    structure_field = c("name")
  )

# exporting
database$writeInterim(data_standard)
