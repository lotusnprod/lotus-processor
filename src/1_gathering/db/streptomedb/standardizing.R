#title: "STREPTOMEDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("streptomedb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(uniqueid,
         name,
         smiles,
         pubchem,
         reference_pubchem = pubmedid,
         biologicalsource) %>%
  cSplit("reference_pubchem", sep = ";", direction = "long") %>%
  cSplit("biologicalsource", sep = ";", direction = "long") %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "str_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_pubchem")
  )

# exporting
database$writeInterim(data_standard)
