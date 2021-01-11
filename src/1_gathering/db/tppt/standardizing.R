# title: "TPPT cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(readxl)
library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("tppt")

## files
data_original_1 <- read_excel(database$sourceFiles$tsv,
  sheet = 1
) %>%
  mutate_all(as.character)

data_original_2 <- read_excel(database$sourceFiles$tsv,
  sheet = 3
) %>%
  mutate_all(as.character)

data_filled <- data_original_1 %>%
  mutate(smiles = ifelse(
    Stereo_SMILES == "NI",
    Canonical_SMILES,
    ifelse(
      Stereo_SMILES == "NS",
      Canonical_SMILES,
      ifelse(Stereo_SMILES == "racemat",
        Canonical_SMILES,
        Stereo_SMILES
      )
    )
  ))

# joining
data_original <- left_join(data_filled, data_original_2)

# selecting
data_selected <- data_original %>%
  select(
    Phytotoxin_number,
    name = Phytotoxin_name,
    CASRN,
    smiles,
    PubChem_CID,
    biologicalsource = Latin_plant_name,
    reference = References
  ) %>%
  cSplit("reference", sep = ",", direction = "long") %>%
  mutate_all(as.character) %>%
  mutate(
    reference_external = ifelse(
      test = reference == "clinitox.ch" |
        reference == "KNApSAcK Database" |
        reference == "KNApSAcKDatabase" |
        reference == "EFSA Reoport (2012)" |
        reference == "EFSA Report" |
        reference == "EFSA Report (2010)" |
        reference == "EFSA Report (2012)",
      yes = reference,
      no = NA
    ),
    reference_authors = gsub(
      pattern = "clinitox.ch",
      replacement = "",
      x = reference,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "KNApSAcK Database",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "KNApSAcKDatabase",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "EFSA Reoport (2012)",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "EFSA Report",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "EFSA Report (2010)",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_authors = gsub(
      pattern = "EFSA Report (2012)",
      replacement = "",
      x = reference_authors,
      fixed = TRUE
    )
  ) %>%
  data.frame()

data_selected$reference_authors <-
  y_as_na(x = data_selected$reference_authors, y = "")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tpp_1",
    structure_field = c("name", "smiles"),
    reference_field = c("reference_authors", "reference_external")
  )

# exporting
database$writeInterim(data_standard)
