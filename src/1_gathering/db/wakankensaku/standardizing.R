# title: "wakankensaku cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("wakankensaku")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
)

# manipulating
data_manipulated <- data_original %>%
  select(
    name = Compound,
    biologicalsource = `Plant resources`,
    reference_original = Literature
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "wak_1",
    structure_field = "name",
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)
