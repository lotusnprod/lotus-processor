# title: "Biophytmol cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("biophytmol")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(uniqueid,
         name,
         smiles,
         biologicalsource,
         reference) %>%
  cSplit("biologicalsource", "     ") %>%
  select(uniqueid,
         name,
         smiles,
         biologicalsource = biologicalsource_1,
         reference) %>%
  mutate_all(as.character) %>%
  tibble()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_2",
    structure_field = c("name", "smiles")
  )

# exporting
database$writeInterim(data_standard)

