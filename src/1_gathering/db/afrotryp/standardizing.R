# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")
source("r/database.R")
source("r/standardizing_original.R")

groundhog.library(Hmisc, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

# get paths
database <- databases$get("afrotryp")

# files
data_original <- vroom(
  file = unz(database$sourceFiles$tsv, "AFROTRYP.tsv"),
  delim = "\t",
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

data_selected$biologicalsource <-
  capitalize(tolower(data_selected$biologicalsource))

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