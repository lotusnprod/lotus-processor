# title: "PROCARDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("procardb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
## atomizing ref
data_selected <- data_original %>%
  mutate(
    reference = gsub("(\\d+\\.)([[:alpha:]])", "| \\2", REFERENCES),
    reference = gsub("(\\d+\\.)(\\s)([[:alpha:]])", "| \\3", reference),
    reference = gsub("(\\d+\\.)(\\s)(\\s)([[:alpha:]])", "| \\4", reference),
    reference = sub("\\| ", "", reference)
  ) %>%
  select(
    uniqueid = column_label,
    name = `CAROTENOID NAME`,
    biologicalsource,
    #inchi = InChI #is an inchikey!!!
    smiles = `CANONICAL SMILES`,
    pubchem = `PUBCHEM ID`,
    reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "pro_1",
    structure_field = c("name", "smiles")
  )

# exporting
database$writeInterim(data_standard)
