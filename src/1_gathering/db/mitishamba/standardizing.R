# title: "MITISHAMBA cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(Hmisc)
library(splitstackshape)
library(tidyverse)

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
    reference_original = authors
  ) %>%
  mutate(reference_test = sub(
    pattern = "\\([0-9]{4}\\)",
    replacement = "§",
    x = reference_original
  )) %>%
  cSplit("reference_test", sep = "§") %>%
  mutate(reference_test_2 = sub(
    pattern = "^\\.",
    replacement = "",
    x = reference_test_2
  )) %>%
  mutate(reference_test_2 = sub(
    pattern = "^\\,",
    replacement = "",
    x = reference_test_2
  )) %>%
  mutate(reference_test_2 = trimws(x = reference_test_2)) %>%
  select(
    name,
    biologicalsource,
    smiles,
    reference_authors = reference_test_1,
    reference_original,
    reference_split = reference_test_2
  ) %>%
  data.frame()

data_corrected <- data_selected %>%
  cSplit("biologicalsource", sep = ",", direction = "long") %>%
  filter(grepl(pattern = "[A-Z]", x = biologicalsource)) %>%
  mutate(biologicalsource = capitalize(tolower(biologicalsource))) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected,
    db = "mit_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_original",
      "reference_authors",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)