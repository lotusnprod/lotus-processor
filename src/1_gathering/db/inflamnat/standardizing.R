# title: "INFLAMNAT cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(readxl, date = groundhog.day)
groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)

# get paths
database <- databases$get("inflamnat")

## files
data_original <-
  read_excel(
    path = database$sourceFiles$tsv,
    sheet = 1
  ) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = Index,
    name = Name,
    smiles = SMILES,
    biologicalsource = Origin,
    pubchem = CID,
    reference = Reference
  )

data_manipulated <- data_selected %>%
  mutate(reference_title = gsub(
    pattern = "\"",
    replacement = "",
    x = str_extract(string = reference, pattern = "\".*\"")
  )) %>%
  mutate(
    reference_authors_1 = ifelse(
      test = !is.na(reference_title),
      yes = str_extract(string = reference, pattern = "^[^\\(]+"),
      no = NA
    ),
    reference_publishingDetails = ifelse(
      test = !is.na(reference_title),
      yes = gsub(
        pattern = "\" ",
        replacement = "",
        x = str_extract(string = reference, pattern = "\" .*")
      ),
      no = NA
    ),
    reference_2 = ifelse(
      test = !is.na(reference_title),
      yes = NA,
      no = reference
    )
  ) %>%
  cSplit("reference_2",
    sep = "\\.,",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    reference_title_2 = ifelse(
      test = !is.na(reference_2_3),
      yes = reference_2_3,
      no = reference_2_2
    ),
    reference_authors_2 = ifelse(
      test = !is.na(reference_2_3),
      yes = paste(reference_2_1, reference_2_2, sep = ".,"),
      no = reference_2_1
    )
  ) %>%
  cSplit("reference_authors_2",
    sep = "et.al.", stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(reference_title_2_2 = gsub("\\([0-9]{4}\\)", "", reference_authors_2_2)) %>%
  mutate(
    reference_authors_3 = reference_authors_2_1,
    reference_original = ifelse(
      test = !is.na(reference_title_2_2),
      yes = reference_title_2_2,
      no = reference_title_2
    )
  ) %>%
  mutate(
    reference_authors = ifelse(
      test = !is.na(reference_authors_2_1),
      yes = reference_authors_2_1,
      no = reference_authors_1
    ),
    reference_original = ifelse(
      test = !is.na(reference_title_2_2),
      yes = reference_title_2_2,
      no = reference_title_2
    ),
    reference_original = gsub(
      pattern = "\"",
      replacement = "",
      x = reference_original
    )
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "inf_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_authors",
      "reference_title",
      "reference_original",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
