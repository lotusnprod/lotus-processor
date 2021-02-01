# title: "DrDuke cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("drduke")

## files
data_common <- vroom(
  file = database$sourceFiles$tsvCommon,
  delim = ",",
  col_types = cols(.default = "c")
) %>%
  select(FNFNUM, CNNAM)

data_farmacy <- vroom(
  file = database$sourceFiles$tsvFarmacy,
  delim = ",",
  col_types = cols(.default = "c")
)

data_fntax <- vroom(
  file = database$sourceFiles$tsvTaxa,
  delim = ",",
  col_types = cols(.default = "c")
) %>%
  select(FNFNUM, TAXON)

data_reference <- vroom(
  file = database$sourceFiles$tsvReference,
  delim = ",",
  col_types = cols(.default = "c")
) %>%
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
    reference_external
  )

data_filtered_2 <- data_selected %>%
  filter(!grepl(pattern = "[0-9]", x = REFERENCE)) %>%
  mutate(
    reference_external = reference_original,
    reference_authors = NA,
    reference_original = NA
  ) %>%
  select(-REFERENCE)

data_filtered <- rbind(data_filtered_1, data_filtered_2) %>%
  mutate(reference_external = ifelse(
    test = is.na(reference_external),
    yes = "DRDUKE",
    no = reference_external
  )) %>%
  select(
    organism_clean = biologicalsource,
    structure_name = name,
    everything()
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_filtered,
    db = "duk_1",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_original",
      "reference_external",
      "reference_authors"
    )
  )

# exporting
database$writeInterim(data_standard)
