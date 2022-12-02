# title: "NPCARE cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")
source("r/y_as_na.R")

library(dplyr)
library(readr)
library(stringr)

# get paths
database <- databases$get("npcare")

## files
data_original <- readr::read_delim(
  file = database$sourceFiles$tsv,
  col_types = cols(.default = "c")
)

# selecting
data_selected <- data_original |>
  dplyr::mutate(reference_pubmed = stringr::str_extract(string = ref_link, pattern = "[0-9]{6,9}")) |>
  dplyr::select(
    originalid = id,
    structure_name = compounds,
    structure_smiles = canonical_smiles,
    biologicalpart = extract,
    pubchem = pid,
    organism_clean = species,
    reference_title = ref,
    reference_pubmed
  )

data_selected[] <-
  lapply(data_selected, function(x) {
    gsub(
      pattern = "\"",
      replacement = " ",
      x = x
    )
  })

data_selected <- data_selected %>%
  mutate_all(
    .tbl = .,
    .funs = trimws
  )

data_selected[] <-
  lapply(data_selected, function(x) {
    y_as_na(x, y = "")
  })

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npcare",
    structure_field = c("structure_name", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c("reference_title", "reference_pubmed")
  )

# exporting
database$writeInterim(data_standard)
