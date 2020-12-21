# title: "wakankensaku cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("wakankensaku")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

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
    structure_field = c("name"),
    reference_field = c("reference_original")
  )

# exporting
database$writeInterim(data_standard)
