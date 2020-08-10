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
    reference_original = LONGREF,
    REFERENCE
  )

data_filtered_1 <- data_selected %>%
  filter(grepl(pattern = "[0-9]", x = REFERENCE)) %>%
  mutate(
    reference_unsplittable = sub(
      pattern = "[0-9]{4}.",
      replacement = "§",
      x = reference_original
    )
  ) %>%
  cSplit("reference_unsplittable", sep = "§") %>%
  mutate_all(as.character) %>%
  mutate(
    reference_authors = ifelse(
      test = !is.na(reference_unsplittable_2),
      yes = reference_unsplittable_1,
      no = NA
    ),
    reference_split = reference_unsplittable_2,
    reference_external = NA
  ) %>%
  select(
    name,
    biologicalsource,
    reference_authors,
    reference_original,
    reference_external,
    reference_split
  )

data_filtered_2 <- data_selected %>%
  filter(!grepl(pattern = "[0-9]", x = REFERENCE)) %>%
  mutate(
    reference_external = reference_original,
    reference_authors = NA,
    reference_original = NA,
    reference_split = NA
  ) %>%
  select(-REFERENCE)

data_filtered <- rbind(data_filtered_1, data_filtered_2) %>%
  mutate(reference_external = ifelse(
    test = is.na(reference_external),
    yes = "DRDUKE",
    no = reference_external
  )) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_filtered,
    db = "duk_1",
    structure_field = "name",
    reference_field = c(
      "reference_original",
      "reference_external",
      "reference_authors",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)
