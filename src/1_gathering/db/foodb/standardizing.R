# title: "FOODB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("foodb")

# files
compounds_flavors <-
  readr::read_delim(file = database$sourceFiles$tsvCompoundsFlavors) |>
  dplyr::mutate_all(as.character)

compounds <-
  readr::read_delim(file = database$sourceFiles$tsvCompounds) |>
  dplyr::mutate_all(as.character)

contents <-
  readr::read_delim(file = database$sourceFiles$tsvContent) |>
  dplyr::mutate_all(as.character)

flavors <-
  readr::read_delim(file = database$sourceFiles$tsvFlavor) |>
  dplyr::mutate_all(as.character)

foods <- readr::read_delim(file = database$sourceFiles$tsvFood) |>
  dplyr::mutate_all(as.character)

references <-
  readr::read_delim(file = database$sourceFiles$tsvReference) |>
  dplyr::mutate_all(as.character)

## Compiling flavors
compiled_flavors <- dplyr::full_join(compounds_flavors,
  flavors,
  by = c("flavor_id" = "id"),
  match = "all"
)

clean_flavors <- compiled_flavors |>
  dplyr::select(
    compound_id,
    flavor_id,
    flavor_citations = citations,
    flavor_name = name,
    flavor_group
  )

# Casting
## contents
compounds_contents <- dplyr::left_join(compounds,
  contents,
  by = c("id" = "source_id"),
  match = "all"
) |>
  dplyr::group_by(id) |>
  dplyr::distinct(orig_food_scientific_name,
    orig_food_part,
    .keep_all = TRUE
  ) |>
  dplyr::ungroup()

## flavors
compounds_flavors <- dplyr::left_join(compounds,
  clean_flavors,
  by = c("id" = "compound_id"),
  match = "all"
)

compounds_flavors_casted <- compounds_flavors |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    flavor_id = paste(flavor_id, collapse = "|"),
    flavor_citations = paste(flavor_citations, collapse = "|"),
    flavor_name = paste(flavor_name, collapse = "|"),
    flavor_group = paste(flavor_group, collapse = "|")
  ) |>
  dplyr::ungroup()

compounds_contents_flavors <- dplyr::left_join(
  compounds_contents,
  compounds_flavors_casted
)

# Minimal output
foodb <- compounds_contents_flavors |>
  tidyr::replace_na(list(
    orig_food_scientific_name = "",
    orig_food_common_name = ""
  )) |>
  dplyr::mutate(
    organism_clean = orig_food_scientific_name,
    organism_dirty = orig_food_common_name,
    reference_external = ifelse(
      test = citation_type == "DATABASE" |
        citation_type == "UNKNOWN" |
        citation_type == "EXPERIMENTAL" |
        citation_type == "PREDICTED" |
        citation == "MANUAL",
      yes = citation,
      no = NA
    ),
    reference_pubmed = ifelse(
      test = citation != "MANUAL" &
        citation_type == "ARTICLE" | citation_type == "TEXTBOOK",
      yes = stringr::str_extract(string = citation, pattern = "[0-9]{6,9}"),
      no = NA
    ),
    reference_original = ifelse(
      test =
        !grepl(pattern = "[0-9]{6,9}", x = citation) &
          citation != "MANUAL" &
          citation_type == "ARTICLE" |
          citation_type == "TEXTBOOK" |
          is.na(citation_type),
      yes = citation,
      no = NA
    )
  ) |>
  dplyr::select(
    uniqueid = public_id,
    structure_name = name,
    organism_clean,
    organism_dirty,
    orig_food_part,
    standard_content,
    reference_external,
    reference_pubmed,
    reference_original,
    structure_smiles = moldb_smiles,
    structure_inchi = moldb_inchi,
    structure_inchikey = moldb_inchikey,
    flavor_name,
    flavor_group
  ) |>
  dplyr::mutate(reference_doi = stringr::str_extract(
    pattern = "10.*",
    string = str_extract(
      pattern = "doi.*",
      string = reference_original
    )
  )) |>
  data.frame()

foodb$organism_clean <- trimws(foodb$organism_clean)
foodb$organism_clean <- y_as_na(foodb$organism_clean, "")

foodb$organism_dirty <- trimws(foodb$organism_dirty)
foodb$organism_dirty <- y_as_na(foodb$organism_dirty, "")

# standardizing
data_standard <- standardizing_original(
  data_selected = foodb,
  db = "foodb",
  structure_field = "structure_smiles",
  organism_field = c("organism_clean", "organism_dirty"),
  reference_field = c(
    "reference_external",
    "reference_pubmed",
    "reference_original",
    "reference_doi"
  )
)

# exporting
database$writeInterim(data_standard)
