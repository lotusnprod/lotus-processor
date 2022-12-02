# title: "PhenolExplorer cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(readr)
library(readxl)
library(tidyr)

# get paths
database <- databases$get("phenolexplorer")

# loading all files
compounds_structures <-
  readr::read_delim(file = database$sourceFiles$tsvCompoundsStructures)

compounds <-
  readr::read_delim(file = database$sourceFiles$tsvCompounds)

foods <- readr::read_delim(file = database$sourceFiles$tsvFoods)

publications <-
  readr::read_delim(file = database$sourceFiles$tsvPublications)

composition <-
  readxl::read_excel(
    path = database$sourceFiles$tsvComposition,
    sheet = 1
  )

### joining
a <- dplyr::full_join(compounds, compounds_structures)
b <- dplyr::full_join(a, composition, by = c("name" = "compound"))
c <- dplyr::full_join(b, foods, by = c("food" = "name"))

# pivoting to join right references
colnames(c)[29] <- "publicationids"
colnames(c)[30] <- "pubmedids"

data_pivoted <- c %>%
  splitstackshape::cSplit("publicationids", ";") %>%
  splitstackshape::cSplit("pubmedids", ";") %>%
  dplyr::group_by(id.x) %>%
  tidyr::pivot_longer(
    37:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  dplyr::ungroup()

# adding references
data_referenced <-
  dplyr::left_join(data_pivoted, publications, by = c("publicationids" = "id")) |>
  dplyr::select(
    uniqueid = id.x,
    name,
    cas = cas_number,
    pubchem = pubchem_compound_id,
    smiles,
    standardcontent = mean,
    biologicalsource = food_source_scientific_name,
    reference_title = title,
    reference_pubmed = pubmedids,
    reference_journal = journal_name,
    reference_authors = authors
  ) |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = " and ",
    fixed = TRUE,
    stripWhite = FALSE,
    direction = "long"
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::select(
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    dplyr::everything()
  ) |>
  data.frame()

data_referenced[] <-
  lapply(data_referenced, function(x) {
    gsub(
      pattern = "NULL",
      replacement = NA,
      x = x
    )
  })

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_referenced,
    db = "phenolexplorer",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_title",
      "reference_pubmed",
      "reference_journal",
      "reference_authors"
    )
  )

# exporting
database$writeInterim(data_standard)
