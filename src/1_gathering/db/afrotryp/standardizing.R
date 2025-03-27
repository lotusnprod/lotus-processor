# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")
source("r/capitalize.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("afrotryp")

# files
data_original <-
  readr::read_delim(
    file = unz(
      description = database$sourceFiles$tsv,
      filename = "AFROTRYP.tsv"
    )
  ) |>
  dplyr::mutate_all(as.character)

# selecting
data_selected <- data_original |>
  dplyr::select(
    uniqueid = `Compound code`,
    structure_name = `Compound name`,
    organism_clean = `Species name`,
    biologicalpart = `plant part`,
    reference_authors = Reference,
    reference_publishingDetails = Publisher
  )

data_selected$organism_clean <-
  capitalize(tolower(data_selected$organism_clean))

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "afrotryp",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c("reference_authors", "reference_publishingDetails")
  )

# exporting
database$writeInterim(data_standard)
