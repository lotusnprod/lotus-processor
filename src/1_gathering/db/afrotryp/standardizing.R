# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")
source("r/database.R")
source("r/standardizing_original.R")

library(Hmisc)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("afrotryp")

# files
data_original <-
  vroom(
    file = unz(database$sourceFiles$tsv, "AFROTRYP.tsv"),
    delim = "\t",
  ) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
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
    db = "afr_1",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c("reference_authors", "reference_publishingDetails")
  )

# exporting
database$writeInterim(data_standard)
