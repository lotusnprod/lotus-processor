# title: "NPASS cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("npass")

## files
data_original_1 <- read_delim(
  file = database$sourceFiles$tsvGeneral
) %>%
  mutate_all(as.character) %>%
  data.frame()

data_original_2 <- read_delim(
  file = database$sourceFiles$tsvProperties
) %>%
  mutate_all(as.character) %>%
  data.frame()

data_original_3 <- read_delim(
  file = database$sourceFiles$tsvSpeciesInfo,
  col_types = cols(.default = "c")
) %>%
  data.frame()

data_original_4 <- read_delim(
  file = database$sourceFiles$tsvSpeciesPair,
  col_types = cols(.default = "c")
) %>%
  data.frame()

# joining
data_original <- left_join(data_original_1, data_original_2)

data_original <- left_join(data_original, data_original_4)

data_original <- left_join(data_original, data_original_3)

# selecting
data_selected <- data_original %>%
  select(
    pubchem = pubchem_cid,
    np_id,
    structure_name = pref_name,
    structure_inchi = standard_inchi,
    standard_inchi_key,
    structure_smiles = canonical_smiles,
    organism_clean = org_name,
    reference = ref_id,
    referenceType = ref_id_type
  )

data_manipulated <- data_selected %>%
  mutate(
    reference_doi = ifelse(
      test = referenceType == "DOI",
      yes = reference,
      no = NA
    ),
    reference_pubmed = ifelse(
      test = referenceType == "PMID",
      yes = reference,
      no = NA
    ),
    reference_external = ifelse(
      test = referenceType == "Database" |
        referenceType == "Patent" |
        referenceType == "Dataset",
      yes = reference,
      no = NA
    ),
    reference_title = ifelse(
      test = referenceType == "Publication" | referenceType == "Book",
      yes = reference,
      no = NA
    ),
  ) %>%
  data.frame()

data_manipulated$structure_name <- y_as_na(data_manipulated$structure_name, "n.a.")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "npass",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_doi",
      "reference_pubmed",
      "reference_external",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
