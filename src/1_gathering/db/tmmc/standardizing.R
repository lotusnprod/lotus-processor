# title: "TMMC cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readxl)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("tmmc")

## files
data_original <- readxl::read_excel(
  path = database$sourceFiles$tsv,
  sheet = 1
) |>
  dplyr::mutate_all(as.character)

data_original_long <- data_original %>%
  splitstackshape::cSplit("CSID", "|") %>%
  tidyr::pivot_longer(17:ncol(.)) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::distinct(SCIENCE, COMPOUND, .keep_all = TRUE) %>%
  dplyr::mutate(
    name = COMPOUND,
    biologicalsource = stringr::str_extract(
      string = SCIENCE,
      pattern = "(?<=\\[).+?(?=\\])"
    ),
    reference_pubmed = stringr::str_extract(
      string = LINK,
      pattern = "[0-9]{6,9}"
    )
  ) %>%
  dplyr::distinct(name, .keep_all = TRUE) %>%
  dplyr::mutate(
    biologicalsource = gsub(
      pattern = "<i>",
      replacement = "",
      x = biologicalsource
    ),
    biologicalsource = gsub(
      pattern = "</i>",
      replacement = "",
      x = biologicalsource
    )
  ) %>%
  data.frame() %>%
  dplyr::select(
    structure_name = name,
    organism_clean = biologicalsource,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original_long,
    db = "tmmc",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = "reference_pubmed"
  )

# exporting
database$writeInterim(data_standard)
