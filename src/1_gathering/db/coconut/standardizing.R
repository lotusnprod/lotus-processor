# title: "COCONUT cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("coconut")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  quote = ""
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
    fixed = TRUE
  ) %>%
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
  mutate(n = nchar(reference_doi)) %>%
  mutate(
    reference_pubmed = ifelse(
      test = n == 8,
      yes = str_extract(string = reference_doi, pattern = "[0-9]{8}"),
      no = NA
    ),
    reference_doi = ifelse(
      test = n >= 9,
      yes = str_extract(string = reference_doi, pattern = "^10.*"),
      no = NA
    ),
    reference_publishingDetails = ifelse(
      test = n >= 20,
      yes = str_extract(string = reference_publishingDetails, pattern = "^[A-Z].*"),
      no = NA
    ),
    reference_authors = ifelse(
      test = n < 20,
      yes = str_extract(string = reference_authors, pattern = "^[A-Z].*"),
      no = NA
    )
  ) %>%
  select(-n)

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

findCapitals_1 <- data_corrected$biologicalsource %>%
  data.frame() %>%
  filter(!is.na(.)) %>%
  mutate(nchar = nchar(.)) %>%
  mutate(first = substr(., start = 1, stop = 1)) %>%
  mutate(last = substr(., start = nchar, stop = nchar)) %>%
  select(
    x = ".",
    nchar,
    first,
    last
  ) %>%
  distinct()

findCapitals_2 <- data_corrected$biologicalsource %>%
  data.frame() %>%
  filter(!is.na(.)) %>%
  mutate(nchar = nchar(.)) %>%
  mutate(first = substr(., start = 1, stop = 1)) %>%
  mutate(last = substr(., start = nchar, stop = nchar)) %>%
  select(
    y = ".",
    nchar,
    first,
    last
  ) %>%
  distinct()

findCapitals_3 <- full_join(findCapitals_1, findCapitals_2) %>%
  filter(tolower(x) == tolower(y)) %>%
  filter(x != y) %>%
  mutate(
    capitals_x = str_count(string = x, pattern = "[A-Z]"),
    capitals_y = str_count(string = y, pattern = "[A-Z]")
  ) %>%
  filter(capitals_x > capitals_y)

data_corrected_capitals <- left_join(data_corrected,
  findCapitals_3,
  by = c("biologicalsource" = "x")
) %>%
  mutate(organism_clean = ifelse(
    test = !is.na(y),
    yes = y,
    no = biologicalsource
  )) %>%
  select(
    -y,
    -nchar,
    -first,
    -last,
    -capitals_x,
    -capitals_y
  ) %>%
  mutate(n = str_count(organism_clean, "\\S+")) %>%
  mutate(organism_dirty = ifelse(test = n == 1,
    yes = organism_clean,
    no = NA
  )) %>%
  mutate(organism_clean = ifelse(test = n > 1,
    yes = organism_clean,
    no = NA
  )) %>%
  select(
    structure_inchi = inchi,
    structure_smiles = smiles,
    structure_name = name,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected_capitals,
    db = "coc_1",
    structure_field = c("structure_inchi", "structure_smiles", "structure_name"),
    organism_field = c("organism_clean", "organism_dirty"),
    reference_field = c(
      "reference_doi",
      "reference_pubmed",
      "reference_authors",
      "reference_publishingDetails",
      "reference_external"
    )
  )

# exporting
database$writeInterim(data_standard)
