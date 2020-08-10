# title: "PhenolExplorer cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("phenolexplorer")

# loading all files
compounds_classification <- read_delim(
  file = database$sourceFiles$tsvCompoundsClassification,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds_structures <- read_delim(
  file = database$sourceFiles$tsvCompoundsStructures,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

compounds <- read_delim(
  file = database$sourceFiles$tsvCompounds,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods_classification <- read_delim(
  file = database$sourceFiles$tsvFoodsClassification,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

foods <- read_delim(
  file = database$sourceFiles$tsvFoods,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites_structures <- read_delim(
  file = database$sourceFiles$tsvMetabolitesStructures,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

metabolites <- read_delim(
  file = database$sourceFiles$tsvMetabolites,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

publications <- read_delim(
  file = database$sourceFiles$tsvPublications,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

composition <-
  read_excel(path = database$sourceFiles$tsvComposition,
             sheet = 1)

### joining
a <- full_join(compounds, compounds_structures)
b <- full_join(a, composition, by = c("name" = "compound"))
c <- full_join(b, foods, by = c("food" = "name"))


# pivoting to join right references
colnames(c)[29] <- "publicationids"
colnames(c)[30] <- "pubmedids"

data_pivoted <- c %>%
  cSplit("publicationids", ";") %>%
  cSplit("pubmedids", ";") %>%
  group_by(id.x) %>%
  pivot_longer(
    37:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  ungroup()

# adding references
data_referenced <-
  left_join(data_pivoted, publications, by = c("publicationids" = "id")) %>%
  select(
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
  ) %>%
  cSplit(
    "biologicalsource",
    sep = " and ",
    fixed = TRUE,
    stripWhite = FALSE,
    direction = "long"
  ) %>%
  mutate_all(as.character) %>%
  data.frame()

data_referenced[] <-
  lapply(data_referenced, function(x)
    gsub("NULL", NA, x))

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_referenced,
    db = "phe_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_title",
      "reference_pubmed",
      "reference_journal",
      "reference_authors"
    )
  )

# exporting
database$writeInterim(data_standard)
