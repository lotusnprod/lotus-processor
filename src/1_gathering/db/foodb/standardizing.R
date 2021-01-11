# title: "FOODB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("foodb")

# files
compounds_flavors <-
  vroom(
    file = database$sourceFiles$tsvCompoundsFlavors,
    delim = ","
  ) %>%
  mutate_all(as.character)

compounds <- read_delim(
  file = database$sourceFiles$tsvCompounds,
  delim = ","
) %>%
  mutate_all(as.character)

contents <- vroom(
  file = database$sourceFiles$tsvContent,
  delim = ","
) %>%
  mutate_all(as.character)

flavors <- vroom(
  file = database$sourceFiles$tsvFlavor,
  delim = ","
) %>%
  mutate_all(as.character)

foods <- vroom(
  file = database$sourceFiles$tsvFood,
  delim = ",",
  quote = ""
) %>%
  mutate_all(as.character)

references <- vroom(
  file = database$sourceFiles$tsvReference,
  delim = ","
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
    biologicalsource = paste(orig_food_scientific_name,
      orig_food_common_name,
      sep = " "
    ),
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
    name,
    biologicalsource,
    orig_food_part,
    standard_content,
    reference_external,
    reference_pubmed,
    reference_original,
    smiles = moldb_smiles,
    inchi = moldb_inchi,
    inchikey = moldb_inchikey,
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

foodb$biologicalsource <- trimws(foodb$biologicalsource)
foodb$biologicalsource <-
  y_as_na(foodb$biologicalsource, "")

# standardizing
data_standard <- standardizing_original(
  data_selected = foodb,
  db = "foo_1",
  structure_field = c("name", "smiles", "inchi"),
  reference_field = c(
    "reference_external",
    "reference_pubmed",
    "reference_original",
    "reference_doi"
  )
)

# exporting
database$writeInterim(data_standard)
