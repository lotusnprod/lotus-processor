# title: "DrDuke cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(readr)

# get paths
database <- databases$get("drduke")

## files
data_common <- readr::read_delim(
  file = database$sourceFiles$tsvCommon,
  delim = ",",
  col_types = cols(.default = "c")
) |>
  dplyr::select(FNFNUM, CNNAM)

data_farmacy <- readr::read_delim(
  file = database$sourceFiles$tsvFarmacy,
  delim = ",",
  col_types = cols(.default = "c")
)

data_fntax <- readr::read_delim(
  file = database$sourceFiles$tsvTaxa,
  delim = ",",
  col_types = cols(.default = "c")
) |>
  dplyr::select(FNFNUM, TAXON)

data_reference <- readr::read_delim(
  file = database$sourceFiles$tsvReference,
  delim = ",",
  col_types = cols(.default = "c")
) |>
  dplyr::select(REFERENCE, LONGREF)

# joining
data_joined <- dplyr::left_join(data_farmacy, data_fntax)

data_joined <- dplyr::left_join(data_joined, data_reference)

# selecting
data_selected <- data_joined |>
  dplyr::select(
    name = CHEM,
    biologicalsource = TAXON,
    reference_original = LONGREF,
    REFERENCE
  )

data_filtered_1 <- data_selected |>
  dplyr::filter(grepl(pattern = "[0-9]", x = REFERENCE)) |>
  dplyr::mutate(
    reference_unsplittable = sub(
      pattern = "[0-9]{4}.",
      replacement = "§",
      x = reference_original
    )
  ) |>
  splitstackshape::cSplit("reference_unsplittable", sep = "§") |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate(
    reference_authors = ifelse(
      test = !is.na(reference_unsplittable_2),
      yes = reference_unsplittable_1,
      no = NA
    ),
    reference_split = reference_unsplittable_2,
    reference_external = NA
  ) |>
  dplyr::select(
    name,
    biologicalsource,
    reference_authors,
    reference_original,
    reference_external
  )

data_filtered_2 <- data_selected |>
  dplyr::filter(!grepl(pattern = "[0-9]", x = REFERENCE)) |>
  dplyr::mutate(
    reference_external = reference_original,
    reference_authors = NA,
    reference_original = NA
  ) |>
  dplyr::select(-REFERENCE)

data_filtered <- rbind(data_filtered_1, data_filtered_2) |>
  dplyr::mutate(reference_external = ifelse(
    test = is.na(reference_external),
    yes = "DRDUKE",
    no = reference_external
  )) |>
  dplyr::select(
    organism_clean = biologicalsource,
    structure_name = name,
    dplyr::everything()
  ) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_filtered,
    db = "drduke",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_original",
      "reference_external",
      "reference_authors"
    )
  )

# exporting
database$writeInterim(data_standard)
