# title: "STREPTOMEDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
ibrary(vroom)

# get paths
database <- databases$get("streptomedb")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  col_names = TRUE,
  id = NULL,
  progress = TRUE,
  escape_double = FALSE,
  trim_ws = TRUE,
  quote = ""
)

## selecting
data_selected <- data_original %>%
  select(uniqueid,
    name,
    smiles,
    pubchem,
    reference_pubmed = pubmedid,
    biologicalsource
  ) %>%
  cSplit("reference_pubmed", sep = ";", direction = "long") %>%
  cSplit("biologicalsource", sep = ";", direction = "long") %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "str_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_pubmed")
  )

# exporting
database$writeInterim(data_standard)
