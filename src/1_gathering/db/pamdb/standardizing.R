# title: "PAMDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readxl)
library(splitstackshape)
library(stringr)

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
    structure_name = Name,
    structure_inchi = InChI,
    structure_smiles = SMILES,
    cas = `CAS number`,
    reference = References
  ) %>%
  mutate(organism_clean = "Pseudomonas aeruginosa")

data_manipulated <- data_selected %>%
  cSplit("reference",
         sep = "Pubmed:",
         fixed = TRUE,
         stripWhite = FALSE) %>%
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
         stripWhite = FALSE) %>%
  mutate_all(as.character) %>%
  select(
    uniqueid,
    organism_clean,
    structure_name,
    structure_inchi,
    structure_smiles,
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
    db = "pamdb",
    structure_field = c("structure_name", "structure_inchi", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c("reference_original",
                        "reference_pubmed",
                        "reference_title")
  )

# exporting
database$writeInterim(data_standard)
