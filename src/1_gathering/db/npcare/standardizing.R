# title: "NPCARE cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")
source("r/y_as_na.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("npcare")

## files
data_original <- vroom(
  file = database$sourceFiles$tsv,
  delim = ";",
  col_types = cols(.default = "c")
)

# selecting
data_selected <- data_original %>%
  mutate(reference_pubmed = str_extract(string = ref_link, pattern = "[0-9]{6,9}")) %>%
  select(
    originalid = id,
    name = compounds,
    smiles = canonical_smiles,
    biologicalpart = extract,
    pubchem = pid,
    biologicalsource = species,
    reference_title = ref,
    reference_pubmed
  )

data_selected[] <-
  lapply(data_selected, function(x) {
    gsub("\"", " ", x)
  })

data_selected <- data_selected %>%
  mutate_all(
    .tbl = .,
    .funs = trimws
  )

data_selected[] <-
  lapply(data_selected, function(x) {
    y_as_na(x, y = "")
  })

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npc_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_title", "reference_pubmed")
  )

# exporting
database$writeInterim(data_standard)
