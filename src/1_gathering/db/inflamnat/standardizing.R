# title: "INFLAMNAT cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readxl)
library(splitstackshape)
library(stringr)

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
    structure_name = Name,
    structure_smiles = SMILES,
    organism_clean = Origin,
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
  mutate(n = str_count(organism_clean, "\\S+")) %>%
  mutate(organism_dirty = ifelse(test = n > 1,
    yes = organism_clean,
    no = NA
  )) %>%
  mutate(organism_clean = ifelse(test = n == 1,
    yes = organism_clean,
    no = NA
  )) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "inflamnat",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = c("organism_clean", "organism_dirty"),
    reference_field = c(
      "reference_authors",
      "reference_title",
      "reference_original",
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
