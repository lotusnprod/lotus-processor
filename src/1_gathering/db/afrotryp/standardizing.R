# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(tidyverse)

# get paths
database <- databases$get("afrotryp")

# files
data_original <- read_delim(
  file = unz(database$sourceFiles$tsv, "AFROTRYP.tsv"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = `Compound code`,
    name = `Compound name`,
    biologicalsource = `Species name`,
    biologicalpart = `plant part`,
    reference_authors = Reference,
    reference_publishingDetails = Publisher
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "afr_1",
    structure_field = "name",
    reference_field = c("reference_authors", "reference_publishingDetails")
  )

# exporting
database$writeInterim(data_standard)