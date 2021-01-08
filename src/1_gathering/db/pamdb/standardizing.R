# title: "PAMDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(readxl, date = groundhog.day)
groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)

# get paths
database <- databases$get("pamdb")

## files
data_original <-
  read_excel(database$sourceFiles$tsv) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = MetID,
    name = Name,
    inchi = InChI,
    smiles = SMILES,
    cas = `CAS number`,
    reference = References
  ) %>%
  mutate(biologicalsource = "Pseudomonas aeruginosa")

data_manipulated <- data_selected %>%
  cSplit("reference",
    sep = "Pubmed:",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    reference_title = str_extract(string = reference_1, pattern = "\".*\""),
    reference_original = ifelse(
      test = !is.na(reference_title),
      yes = NA,
      no = reference_1
    )
  ) %>%
  cSplit("reference_2",
    sep = " ",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  select(
    uniqueid,
    biologicalsource,
    name,
    inchi,
    smiles,
    cas,
    reference_original,
    reference_title,
    reference_pubmed = reference_2_02
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "pam_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c(
      "reference_original",
      "reference_pubmed",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
