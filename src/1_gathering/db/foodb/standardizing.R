# title: "FOODB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("foodb")

# files
compounds_flavors <- read_delim(
  file = database$sourceFiles$tsvCompoundsFlavors,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

compounds <- read_delim(
  file = database$sourceFiles$tsvCompounds,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

contents <- read_delim(
  file = database$sourceFiles$tsvContent,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

flavors <- read_delim(
  file = database$sourceFiles$tsvFlavor,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

foods <- read_delim(
  file = database$sourceFiles$tsvFood,
  delim = ",",
  escape_double = TRUE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

references <- read_delim(
  file = database$sourceFiles$tsvReference,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## Compiling flavors
compiled_flavors <- full_join(compounds_flavors,
                              flavors,
                              by = c('flavor_id' = 'id'),
                              match = "all")

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
                                by = c('id' = 'source_id'),
                                match = "all") %>%
  group_by(id) %>%
  distinct(orig_food_scientific_name,
           orig_food_part,
           .keep_all = TRUE) %>%
  ungroup()

## flavors
compounds_flavors <- left_join(compounds,
                               clean_flavors,
                               by = c('id' = 'compound_id'),
                               match = "all")

compounds_flavors_casted <- compounds_flavors %>%
  group_by(id) %>%
  summarise(
    flavor_id = paste(flavor_id, collapse = "|"),
    flavor_citations = paste(flavor_citations, collapse = "|"),
    flavor_name = paste(flavor_name, collapse = "|"),
    flavor_group = paste(flavor_group, collapse = "|")
  ) %>%
  ungroup()

compounds_contents_flavors <- left_join(compounds_contents,
                                        compounds_flavors_casted)

# Minimal output
foodb <- compounds_contents_flavors %>%
  replace_na(list(
    orig_food_scientific_name = "",
    orig_food_common_name = ""
  )) %>%
  mutate(
    biologicalsource = paste(orig_food_scientific_name,
                             orig_food_common_name,
                             sep = " "),
    reference = paste(citation, citation_type, sep = "§")
  ) %>%
  select(
    uniqueid = public_id,
    name,
    biologicalsource,
    orig_food_part,
    standard_content,
    reference,
    smiles = moldb_smiles,
    inchi = moldb_inchi,
    inchikey = moldb_inchikey,
    flavor_name,
    flavor_group
  )

foodb$biologicalsource <- trimws(foodb$biologicalsource)
foodb$biologicalsource <-
  y_as_na(foodb$biologicalsource, "")
foodb$reference <-
  y_as_na(foodb$reference, "NA NA")

# standardizing
data_standard <- standardizing_original(
  data_selected = foodb,
  db = "foo_1",
  structure_field = c("name", "smiles", "inchi")
)

# exporting
database$writeInterim(data_standard)