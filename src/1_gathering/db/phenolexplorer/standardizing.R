# title: "PhenolExplorer cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(readxl)
library(vroom)

# get paths
database <- databases$get("phenolexplorer")

# loading all files
compounds_classification <-
  vroom(
    file = database$sourceFiles$tsvCompoundsClassification,
    delim = ","
  )

compounds_structures <-
  vroom(
    file = database$sourceFiles$tsvCompoundsStructures,
    delim = ","
  )

compounds <- vroom(
  file = database$sourceFiles$tsvCompounds,
  delim = ","
)

foods_classification <-
  vroom(
    file = database$sourceFiles$tsvFoodsClassification,
    delim = ","
  )

foods <- vroom(
  file = database$sourceFiles$tsvFoods,
  delim = ","
)

metabolites_structures <-
  vroom(
    file = database$sourceFiles$tsvMetabolitesStructures,
    delim = ","
  )

metabolites <- vroom(
  file = database$sourceFiles$tsvMetabolites,
  delim = ","
)

publications <- vroom(
  file = database$sourceFiles$tsvPublications,
  delim = ","
)

composition <-
  read_excel(
    path = database$sourceFiles$tsvComposition,
    sheet = 1
  )

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
  select(
    structure_name = name,
    structure_smiles = smiles,
    organism_clean = biologicalsource,
    everything()
  ) %>%
  data.frame()

data_referenced[] <-
  lapply(data_referenced, function(x) {
    gsub("NULL", NA, x)
  })

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_referenced,
    db = "phe_1",
    structure_field = c("structure_name", "structure_smiles"),
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
