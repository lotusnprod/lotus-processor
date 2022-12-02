# title:"NPEDIA cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("npedia")

## files
data_original <- readr::read_delim(
  file = gzfile(description = database$sourceFiles$tsv),
  col_types = cols(.default = "c")
)

# selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid = ID,
    structure_name = Name,
    organism_clean = Source,
    structure_inchi = InChI,
    cas = `CAS No.`,
    reference_original = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npedia",
    structure_field = c("structure_name", "structure_inchi"),
    organism_field = "organism_clean",
    reference_field = "reference_original"
  )

# exporting
database$writeInterim(data_standard)
