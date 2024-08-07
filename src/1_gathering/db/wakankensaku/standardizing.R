# title: "wakankensaku cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("wakankensaku")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv))

# manipulating
data_manipulated <- data_original |>
  dplyr::select(
    structure_name = Compound,
    organism_dirty = `Plant resources`,
    reference_original = Literature
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "wakankensaku",
    structure_field = "structure_name",
    organism_field = "organism_dirty",
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)
