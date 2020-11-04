# title: "NPCARE cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("npcare")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

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