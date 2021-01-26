# title: "NPASS cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("npass")

## files
data_original_1 <- vroom(
  file = database$sourceFiles$tsvGeneral,
  delim = "\t",
  quote = ""
) %>%
  mutate_all(as.character) %>%
  data.frame()

data_original_2 <- vroom(
  file = database$sourceFiles$tsvProperties,
  delim = "\t"
) %>%
  mutate_all(as.character) %>%
  data.frame()

data_original_3 <- vroom(
  file = database$sourceFiles$tsvSpeciesInfo,
  delim = "\t",
  col_types = cols(.default = "c")
) %>%
  data.frame()

data_original_4 <- vroom(
  file = database$sourceFiles$tsvSpeciesPair,
  delim = "\t",
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
    name = pref_name,
    inchi = standard_inchi,
    standard_inchi_key,
    smiles = canonical_smiles,
    biologicalsource = org_name,
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

data_manipulated$name <- y_as_na(data_manipulated$name, "n.a.")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "npa_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c(
      "reference_doi",
      "reference_pubmed",
      "reference_external",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
