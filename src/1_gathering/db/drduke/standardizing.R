# title: "DrDuke cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("drduke")

## files
data_common <- read_delim(
  file = database$sourceFiles$tsvCommon,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(FNFNUM, CNNAM)

data_farmacy <- read_delim(
  file = database$sourceFiles$tsvFarmacy,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_fntax <- read_delim(
  file = database$sourceFiles$tsvTaxa,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(FNFNUM, TAXON)

data_reference <- read_delim(
  file = database$sourceFiles$tsvReference,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(REFERENCE, LONGREF)

# joining
data_joined <- left_join(data_farmacy, data_fntax)

data_joined <- left_join(data_joined, data_reference)

# selecting
data_selected <- data_joined %>%
  select(
    name = CHEM,
    biologicalsource = TAXON,
    reference_unsplittable = LONGREF
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "duk_1",
    structure_field = "name",
    reference_field = "reference_unsplittable"
  )

# exporting
database$writeInterim(data_standard)
