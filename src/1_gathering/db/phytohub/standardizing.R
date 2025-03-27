# title: "Phytohub cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("phytohub")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv)) |>
  dplyr::mutate_all(as.character)

data_manipulated <- data_original |>
  dplyr::filter(!is.na(biologicalsource)) |>
  dplyr::mutate(
    reference = gsub(
      pattern = "doi: ",
      replacement = "$",
      x = reference
    )
  ) |>
  splitstackshape::cSplit(
    "reference",
    sep = "$",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::filter(grepl(pattern = "^10.", x = reference)) |>
  splitstackshape::cSplit(
    "reference",
    sep = " ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::filter(
    grepl(pattern = "^10.", x = reference) |
      grepl(pattern = "[PubMed:", x = reference, fixed = TRUE)
  ) |>
  dplyr::mutate(
    reference_pubmed = gsub(
      pattern = "[PubMed:",
      replacement = "",
      x = reference,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_pubmed = gsub(
      pattern = "]",
      replacement = "",
      x = reference_pubmed,
      fixed = TRUE
    ),
    reference_doi = reference
  ) |>
  dplyr::select(-reference) |>
  dplyr::mutate(
    reference_doi = ifelse(
      test = grepl(pattern = "^10.", x = reference_doi),
      yes = gsub(
        pattern = "\\.$",
        replacement = "",
        x = reference_doi
      ),
      no = NA_character_
    ),
    reference_pubmed = ifelse(
      test = grepl(pattern = "^10.", x = reference_doi),
      no = reference_pubmed,
      yes = NA_character_
    )
  ) |>
  dplyr::select(
    structure_smiles = smiles,
    structure_inchi = inchi,
    organism_dirty = biologicalsource,
    reference_pubmed,
    reference_doi
  ) |>
  tidyr::pivot_longer(cols = 4:5) |>
  dplyr::filter(!is.na(value)) |>
  tidyr::pivot_wider(names_from = name, values_from = value) |>
  tidyr::unnest_longer(reference_doi) |>
  tidyr::unnest_longer(reference_pubmed)


# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "phytohub",
    structure_field = "structure_inchi",
    organism_field = "organism_dirty",
    reference_field = c(
      "reference_doi",
      "reference_pubmed"
    )
  )

# exporting
database$writeInterim(data_standard)
