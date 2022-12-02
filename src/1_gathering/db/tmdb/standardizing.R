# title: "TMDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("tmdb")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv))

# pivoting
data_pivoted <- data_original %>%
  dplyr::mutate(level = as.numeric(gl(nrow(.) / 28, 28))) %>%
  dplyr::group_by(level) %>%
  tidyr::pivot_wider(names_from = 1, values_from = 2) %>%
  tidyr::unnest() %>%
  dplyr::ungroup()

# selecting
data_selected <- data_pivoted |>
  dplyr::select(
    name = `Entry name`,
    biologicalsource = `Latin name`,
    reference_publishingDetails = References
  ) |>
  splitstackshape::cSplit("reference_publishingDetails",
    sep = ";",
    direction = "long"
  ) |>
  dplyr::mutate(
    reference_publishingDetails = gsub(
      pattern = "\\(.*\\D.*\\)",
      replacement = "",
      x = reference_publishingDetails
    )
  ) |>
  dplyr::select(
    structure_name = name,
    organism_clean = biologicalsource,
    everything()
  ) |>
  data.frame()

data_selected[] <-
  lapply(data_selected, function(x) {
    gsub(
      pattern = "Not Available",
      replacement = NA,
      x = x
    )
  })

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tmdb",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = "reference_publishingDetails"
  )

# exporting
database$writeInterim(data_standard)
