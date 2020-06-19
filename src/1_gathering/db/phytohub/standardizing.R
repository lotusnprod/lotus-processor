# title: "Phytohub cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("phytohub")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = TRUE,
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original,
    db = "phy_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
database$writeInterim(data_standard)
