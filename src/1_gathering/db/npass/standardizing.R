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
data_original_1 <-
  readr::read_delim(file = database$sourceFiles$tsvGeneral) |>
  dplyr::mutate_all(as.character) |>
  data.frame()

data_original_2 <-
  readr::read_delim(
    file = database$sourceFiles$tsvStructure,
    col_names = FALSE
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::select(
    np_id = X1,
    standard_inchi = X2,
    standard_inchi_key = X3,
    canonical_smiles = X4
  ) |>
  data.frame()

data_original_3 <-
  readr::read_delim(
    file = database$sourceFiles$tsvSpeciesInfo,
    col_types = cols(.default = "c")
  ) |>
  data.frame()

data_original_4 <-
  readr::read_delim(
    file = database$sourceFiles$tsvSpeciesPair,
    col_types = cols(.default = "c")
  ) |>
  data.frame()

# joining
data_original <- dplyr::left_join(data_original_1, data_original_2)

data_original <- dplyr::left_join(data_original, data_original_4)

data_original <- dplyr::left_join(data_original, data_original_3)

# selecting
data_selected <- data_original |>
  dplyr::select(
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

data_manipulated <- data_selected |>
  dplyr::mutate(
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
  ) |>
  data.frame()

data_manipulated$structure_name <-
  y_as_na(data_manipulated$structure_name, "n.a.")

# removing completely aberrant structure failing with RDKit
data_manipulated <- data_manipulated |>
  dplyr::filter(
    !grepl(
      pattern = "c1cc(ccc1[CH2+]1C(=Cc2c(cc(cc2O1)O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@@H](CO)O1)O)O)O)O",
      x = structure_smiles,
      fixed = TRUE
    )
  )

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
