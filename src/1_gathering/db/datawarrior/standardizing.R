# title: "datawarrior cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)

# get paths
database <- databases$get("datawarrior")

## files
data_original <-
  readr::read_delim(file = database$sourceFiles$tsv) |>
  dplyr::mutate_all(as.character)

# manipulating
data_manipulated <- data_original |>
  splitstackshape::cSplit("reference", sep = "; ") |>
  splitstackshape::cSplit("remarks",
    sep = "; ",
    direction = "long",
    fixed = TRUE
  ) |>
  splitstackshape::cSplit(
    "remarks",
    sep = "<NL>",
    stripWhite = FALSE,
    direction = "wide",
    fixed = TRUE
  ) |>
  dplyr::select(
    structure_name = name,
    structure_smiles = Smiles,
    structure_inchi = InChI,
    biologicalsource = remarks_1,
    reference_authors = reference_1,
    reference_title = reference_2
  ) |>
  dplyr::mutate(
    organism_clean = gsub(
      pattern = "Origin:",
      replacement = "",
      x = biologicalsource,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate_all(trimws) |>
  dplyr::mutate(
    reference_title = ifelse(
      test = reference_title == "Giftpflanzen Pflanzengifte4. ?berarbeitete Aufl.",
      yes = "Giftpflanzen Pflanzengifte, 4. ?berarbeitete Aufl.",
      no = reference_title
    )
  ) |>
  dplyr::mutate(
    reference_isbn = ifelse(
      test = reference_title == "Giftpflanzen Pflanzengifte, 4. ?berarbeitete Aufl.",
      yes = "978-3-933203-31-1",
      no = "3-8047-1053-0"
    )
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "datawarrior",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c("reference_authors", "reference_title", "reference_isbn")
  )

# exporting
database$writeInterim(data_standard)
