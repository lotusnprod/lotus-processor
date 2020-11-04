# title: "PROCARDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("procardb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
## atomizing ref
data_selected <- data_original %>%
  mutate(
    reference = gsub("(\\d+\\.)([[:alpha:]])", "| \\2", REFERENCES),
    reference = gsub("(\\d+\\.)(\\s)([[:alpha:]])", "| \\3", reference),
    reference = gsub("(\\d+\\.)(\\s)(\\s)([[:alpha:]])", "| \\4", reference),
    reference = sub("\\| ", "", reference)
  ) %>%
  select(
    uniqueid = column_label,
    name = `CAROTENOID NAME`,
    biologicalsource,
    # inchi = InChI #is an inchikey!!!
    smiles = `CANONICAL SMILES`,
    pubchem = `PUBCHEM ID`,
    reference
  )

data_manipulated <- data_selected %>%
  cSplit(
    "reference",
    sep = "|",
    fixed = TRUE,
    stripWhite = FALSE,
    direction = "long",
    drop = FALSE
  ) %>%
  mutate(reference = sub(pattern = ":", replacement = "§", reference)) %>%
  cSplit("reference", sep = "§") %>%
  cSplit("reference_2",
    sep = "PMID",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  mutate_all(as.character) %>%
  mutate(
    reference_authors = reference_1,
    reference_split = reference_2_1,
    reference_pubmed = str_extract(string = reference_2_2, pattern = "[0-9]{6,9}")
  ) %>%
  mutate(reference_split = ifelse(
    test = !grepl(pattern = "[A-Za-z]", x = reference_split),
    yes = reference_authors,
    no = reference_split
  )) %>%
  mutate(
    reference_authors = ifelse(
      test = reference_authors == reference_split,
      yes = NA,
      no = reference_authors
    )
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "pro_1",
    structure_field = c("name", "smiles"),
    reference_field = c(
      "reference_authors",
      "reference_split",
      "reference_pubmed"
    )
  )

# exporting
database$writeInterim(data_standard)