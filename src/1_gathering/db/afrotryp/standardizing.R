# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)

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
  mutate(reference = paste(Reference, Publisher, sep = " ")) %>%
  select(
    uniqueid = `Compound code`,
    name = `Compound name`,
    biologicalsource = `Species name`,
    biologicalpart = `plant part`,
    reference
  )

data_selected$reference <- gsub("NA", "", data_selected$reference)

# standardizing
data_standard <-
  standardizing_original(data_selected = data_selected,
                         db = "afr_1",
                         structure_field = "name")

# exporting
database$writeInterim(data_standard)