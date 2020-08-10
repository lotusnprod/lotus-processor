# title:"NPEDIA cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("npedia")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = TRUE,
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

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
    reference_field = c("reference_original")
  )

# exporting
database$writeInterim(data_standard)
