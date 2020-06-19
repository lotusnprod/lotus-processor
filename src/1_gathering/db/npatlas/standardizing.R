# title: "NPATLAS cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("npatlas")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  mutate(organism = paste(Genus,
                          `Origin Species`,
                          sep = " "),
         reference = `Isolation Reference DOI`) %>%
  select(
    NPAID,
    name = Names,
    InChIKey,
    inchi = InChI,
    smiles = SMILES,
    biologicalsource = organism,
    reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_2",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
database$writeInterim(data_standard)
