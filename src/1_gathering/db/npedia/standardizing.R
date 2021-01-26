# title:"NPEDIA cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("npedia")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  quote = "",
  col_types = cols(.default = "c")
)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    name = Name,
    biologicalsource = Source,
    inchi = InChI,
    cas = `CAS No.`,
    reference_original = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npe_1",
    structure_field = c("name", "inchi"),
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)