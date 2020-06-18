# title: "COCONUT cleaneR"

# loading paths
source("paths.R")
source("functions.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("coconut")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
data_selected <- data_original %>%
  select(
    name,
    inchi = inchi,
    smiles = SMILES,
    biologicalsource = textTaxa,
    reference = citationDOI,
    externalDB = found_in_databases
  ) %>%
  mutate(
    biologicalsource = gsub("\\[", "", biologicalsource),
    biologicalsource = gsub("\\]", "", biologicalsource),
    biologicalsource = gsub("notax", NA, biologicalsource),
    biologicalsource = gsub("\"", "", biologicalsource),
    reference = gsub("\\[", "", reference),
    reference = gsub("\\]", "", reference),
    reference = gsub(",", "|", reference),
    reference = gsub("\"", "", reference),
    externalDB = gsub("\\[", "", externalDB),
    externalDB = gsub("\\]", "", externalDB),
    externalDB = gsub(",", "|", externalDB),
    externalDB = gsub("\"", "", externalDB),
  )

data_selected$name <- y_as_na(data_selected$name, "")
data_selected$inchi <- y_as_na(data_selected$inchi, "")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "")
data_selected$reference <- y_as_na(data_selected$reference, "")
data_selected$externalDB <- y_as_na(data_selected$externalDB, "")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "coc_1",
    structure_field = c("inchi", "smiles", "name")
  )

# exporting
database$writeInterim(data_standard)
