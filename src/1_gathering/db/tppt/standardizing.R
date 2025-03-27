# title: "TPPT cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readxl)
library(splitstackshape)

# get paths
database <- databases$get("tppt")

## files
data_original_1 <-
  readxl::read_excel(
    path = database$sourceFiles$tsv,
    sheet = 1
  ) |>
  dplyr::mutate_all(as.character)

data_original_2 <- readxl::read_excel(database$sourceFiles$tsv, sheet = 3) |>
  dplyr::mutate_all(as.character)

data_filled <- data_original_1 |>
  dplyr::mutate(
    smiles = ifelse(
      Stereo_SMILES == "NI",
      Canonical_SMILES,
      ifelse(
        Stereo_SMILES == "NS",
        Canonical_SMILES,
        ifelse(Stereo_SMILES == "racemat", Canonical_SMILES, Stereo_SMILES)
      )
    )
  )

# joining
data_original <- dplyr::left_join(data_filled, data_original_2)

# selecting
data_selected <- data_original |>
  dplyr::select(
    Phytotoxin_number,
    name = Phytotoxin_name,
    CASRN,
    smiles,
    PubChem_CID,
    biologicalsource = Latin_plant_name,
    reference = References
  ) |>
  splitstackshape::cSplit("reference", sep = ",", direction = "long") |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate(
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
    reference_original = gsub(
      pattern = "clinitox.ch",
      replacement = "",
      x = reference,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "KNApSAcK Database",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "KNApSAcKDatabase",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "EFSA Reoport (2012)",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "EFSA Report",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "EFSA Report (2010)",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_original = gsub(
      pattern = "EFSA Report (2012)",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) |>
  data.frame() |>
  dplyr::select(
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    dplyr::everything()
  )

data_selected$reference_original <-
  y_as_na(x = data_selected$reference_original, y = "")

data_selected$reference_authors <- gsub(
  pattern = "\\([0-9]{4}\\)",
  replacement = "",
  x = data_selected$reference_original
)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tppt",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_original",
      "reference_authors",
      "reference_external"
    )
  )

# exporting
database$writeInterim(data_standard)
