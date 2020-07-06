# title: "MITISHAMBA cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(readxl)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("mitishamba")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  trim_ws = TRUE,
  escape_backslash = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    smiles,
    biologicalsource = plant_species,
    name = common_name,
    reference_unsplittable = authors
  ) %>% 
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "mit_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_unsplittable")
  )

# exporting
database$writeInterim(data_standard)
