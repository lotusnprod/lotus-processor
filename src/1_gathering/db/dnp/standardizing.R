# title: "DNP cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("dnp")

## files
data_original <- vroom(
  file = database$sourceFiles$tsv,
  delim = ","
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = CRC_Number,
    name = Molecule_Name,
    inchi = MolfileName,
    biologicalsource = Biological_Source
  ) %>%
  distinct(uniqueid, .keep_all = TRUE) %>%
  mutate(reference_external = "DNP") %>%
  data.frame()

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "dnp_1",
    structure_field = c("name", "inchi"),
    reference_field = c("reference_external")
  )

# exporting
database$writeInterim(data_standard)
