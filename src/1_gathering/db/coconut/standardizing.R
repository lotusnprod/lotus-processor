# title: "COCONUT cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("coconut")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
data_selected <- data_original %>%
  select(
    name,
    inchi = inchi,
    smiles = SMILES,
    biologicalsource = textTaxa,
    reference_doi = citationDOI,
    reference_external = found_in_databases
  ) %>%
  cSplit(
    "biologicalsource",
    sep = ",",
    direction = "long",
    fixed = TRUE
  ) %>%
  cSplit("reference_doi",
         sep = ",",
         direction = "long",
         fixed = TRUE) %>%
  # cSplit("reference_external",
  #        sep = ",",
  #        direction = "long",
  #        fixed = TRUE) %>%
  mutate(
    biologicalsource = gsub("\\[", "", biologicalsource),
    biologicalsource = gsub("\\]", "", biologicalsource),
    biologicalsource = gsub("notax", "", biologicalsource),
    biologicalsource = gsub("\"", "", biologicalsource),
    reference_doi = gsub("\\[", "", reference_doi),
    reference_doi = gsub("\\]", "", reference_doi),
    reference_doi = gsub("\"", "", reference_doi),
    # reference_external = gsub("\\[", "", reference_external),
    # reference_external = gsub("\\]", "", reference_external),
    # reference_external = gsub("\"", "", reference_external)
  ) %>%
  data.frame()

data_corrected <- data_selected %>%
  mutate(
    reference_publishingDetails = reference_doi,
    reference_authors = reference_doi
  ) %>%
  mutate(
    reference_doi = str_extract(string = reference_doi, pattern = "^10.*"),
    reference_publishingDetails = ifelse(
      test = str_count(string = reference_publishingDetails) >= 20,
      yes = str_extract(string = reference_publishingDetails, pattern = "^[A-Z].*"),
      no = NA
    ),
    reference_authors = ifelse(
      test = str_count(string = reference_authors) < 20,
      yes = str_extract(string = reference_authors, pattern = "^[A-Z].*"),
      no = NA
    )
  )

data_corrected$name <- y_as_na(data_corrected$name, "")
data_corrected$inchi <- y_as_na(data_corrected$inchi, "")
data_corrected$biologicalsource <-
  y_as_na(data_corrected$biologicalsource, "")
data_corrected$reference_external <-
  y_as_na(data_corrected$reference_external, "")

data_corrected$name <- y_as_na(data_corrected$name, "NA")
data_corrected$inchi <- y_as_na(data_corrected$inchi, "NA")
data_corrected$biologicalsource <-
  y_as_na(data_corrected$biologicalsource, "NA")
data_corrected$reference_authors <-
  y_as_na(data_corrected$reference_authors, "NA")
data_corrected$reference_external <-
  y_as_na(data_corrected$reference_external, "NA")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected,
    db = "coc_1",
    structure_field = c("inchi", "smiles", "name"),
    reference_field = c("reference_doi", "reference_authors", "reference_publishingDetails" ,"reference_external")
  )

# exporting
database$writeInterim(data_standard)
