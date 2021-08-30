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
  read_delim(
    file = database$sourceFiles$tsvCompoundsFlavors
  ) %>%
  mutate_all(as.character)

compounds <- read_delim(
  file = database$sourceFiles$tsvCompounds
) %>%
  mutate_all(as.character)

contents <- read_delim(
  file = database$sourceFiles$tsvContent
) %>%
  mutate_all(as.character)

flavors <- read_delim(
  file = database$sourceFiles$tsvFlavor
) %>%
  mutate_all(as.character)

foods <- read_delim(
  file = database$sourceFiles$tsvFood
) %>%
  mutate_all(as.character)

references <- read_delim(
  file = database$sourceFiles$tsvReference,
) %>%
  mutate_all(as.character)

## Compiling flavors
compiled_flavors <- full_join(compounds_flavors,
  flavors,
  by = c("flavor_id" = "id"),
  match = "all"
)

clean_flavors <- compiled_flavors %>%
  select(
    compound_id,
    flavor_id,
    flavor_citations = citations,
    flavor_name = name,
    flavor_group
  )

# Casting
## contents
compounds_contents <- left_join(compounds,
  contents,
  by = c("id" = "source_id"),
  match = "all"
) %>%
  group_by(id) %>%
  distinct(orig_food_scientific_name,
    orig_food_part,
    .keep_all = TRUE
  ) %>%
  ungroup()

## flavors
compounds_flavors <- left_join(compounds,
  clean_flavors,
  by = c("id" = "compound_id"),
  match = "all"
)

compounds_flavors_casted <- compounds_flavors %>%
  group_by(id) %>%
  summarise(
    flavor_id = paste(flavor_id, collapse = "|"),
    flavor_citations = paste(flavor_citations, collapse = "|"),
    flavor_name = paste(flavor_name, collapse = "|"),
    flavor_group = paste(flavor_group, collapse = "|")
  ) %>%
  ungroup()

compounds_contents_flavors <- left_join(
  compounds_contents,
  compounds_flavors_casted
)

# Minimal output
foodb <- compounds_contents_flavors %>%
  replace_na(list(
    orig_food_scientific_name = "",
    orig_food_common_name = ""
  )) %>%
  mutate(
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
      yes = str_extract(string = citation, pattern = "[0-9]{6,9}"),
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
  ) %>%
  select(
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
  ) %>%
  mutate(reference_doi = str_extract(
    pattern = "10.*",
    string = str_extract(
      pattern = "doi.*",
      string = reference_original
    )
  )) %>%
  data.frame()

foodb$organism_clean <- trimws(foodb$organism_clean)
foodb$organism_clean <- y_as_na(foodb$organism_clean, "")

foodb$organism_dirty <- trimws(foodb$organism_dirty)
foodb$organism_dirty <- y_as_na(foodb$organism_dirty, "")

# standardizing
data_standard <- standardizing_original(
  data_selected = foodb,
  db = "foodb",
  structure_field = c("structure_name", "structure_smiles", "structure_inchi"),
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
