# title: "COCONUT cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("coconut")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv)
)

"%ni%" <- Negate("%in%")

# selecting
data_selected <- data_original %>%
  select(
    name,
    inchi = inchi,
    smiles = SMILES,
    biologicalsource = textTaxa,
    reference = citationDOI,
    reference_external = found_in_databases
  ) %>%
  filter(biologicalsource != "[notax]") %>%
  filter(biologicalsource != "[]") %>%
  filter(reference != "[]") %>%
  mutate(
    biologicalsource = gsub("\\[", "", biologicalsource),
    biologicalsource = gsub("\\]", "", biologicalsource),
    biologicalsource = gsub("\"", "", biologicalsource),
    reference = gsub("\\[", "", reference),
    reference = gsub("\\]", "", reference),
    reference = gsub("\"", "", reference),
    # reference_external = gsub("\\[", "", reference_external),
    # reference_external = gsub("\\]", "", reference_external),
    # reference_external = gsub("\"", "", reference_external)
  )

data_counted <- data_selected %>%
  rowwise() %>%
  mutate(
    nrefs = str_count(string = reference, pattern = ", "),
    norganisms = str_count(string = biologicalsource, pattern = ", ")
  )

## right organism, right ref
data_nosplit <- data_counted %>%
  filter(nrefs == 0)

## right organism, right ref
data_easy <- data_counted %>%
  filter(nrefs != 0 & nrefs == norganisms)

## combination of all organisms vs all references since we do not know which one is good
data_hard <- data_counted %>%
  filter(nrefs != 0 & nrefs != norganisms)

data_hard_harborne <- data_hard %>%
  filter(
    grepl(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      x = reference
    )
  )

data_hard_harborne_treated <- data_hard_harborne %>%
  mutate(
    reference = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      replacement = "",
      x = reference
    )
  ) %>%
  cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  group_by(smiles) %>%
  add_count() %>%
  ungroup() %>%
  filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) %>%
  select(-n) %>%
  cSplit("reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  filter(!grepl(pattern = "^,", x = reference)) %>%
  filter(!grepl(pattern = "^\\.", x = reference)) %>%
  filter(!grepl(pattern = "^[0-9]{1-3},", x = reference)) %>%
  filter(!grepl(pattern = "^\\([0-9]{4}\\)", x = reference)) %>%
  filter(str_length(reference) >= 8) %>%
  add_count(reference) %>%
  filter(
    n <= 5400 |
      grepl(
        pattern = "Molecules_2015;20(8):15330-42",
        x = reference
      ) |
      grepl(
        pattern = "Lu,Phytochem.,59,(2002),117",
        x = reference
      )
  ) %>%
  select(-nrefs, -norganisms, -n) %>%
  data.frame()

data_hard_2 <- data_hard %>%
  filter(
    !grepl(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      x = reference
    )
  )

data_nosplit_treated <- data_nosplit %>%
  cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  group_by(smiles) %>%
  add_count() %>%
  ungroup() %>%
  filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) %>%
  select(-n) %>%
  cSplit("reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  filter(str_length(reference) >= 8) %>%
  select(-nrefs, -norganisms) %>%
  data.frame()

data_easy_treated <- data_easy %>%
  cSplit("biologicalsource",
    sep = ", ",
    fixed = TRUE
  ) %>%
  cSplit("reference",
    sep = ", ",
    fixed = TRUE
  ) %>%
  pivot_longer(
    cols = 7:(ncol(.)),
    names_to = c("type", "number"),
    names_sep = "_"
  ) %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  group_by(inchi) %>%
  fill(reference, .direction = "downup") %>%
  group_by(inchi) %>%
  add_count() %>%
  ungroup() %>%
  filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) %>%
  select(-n) %>%
  filter(str_length(reference) >= 8) %>%
  select(-nrefs, -norganisms, -number) %>%
  data.frame()

data_hard_treated <- data_hard_2 %>%
  cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  group_by(smiles) %>%
  add_count() %>%
  ungroup() %>%
  filter(n == 1 | biologicalsource %ni% c("plants", "marine", "fungi")) %>%
  select(-n) %>%
  cSplit(
    "reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  filter(str_length(reference) >= 8) %>%
  add_count(reference) %>%
  filter(n <= 15000) %>%
  select(-nrefs, -norganisms, -n) %>%
  data.frame()

data_treated <- rbind(
  data_nosplit_treated,
  data_easy_treated,
  data_hard_harborne_treated,
  data_hard_treated
) %>%
  filter(grepl(pattern = "[0-9]", x = reference))

data_corrected <- data_treated %>%
  mutate(
    reference_split = reference,
    reference_authors = reference
  ) %>%
  mutate(n = nchar(reference)) %>%
  mutate(
    reference_pubmed = ifelse(
      test = n == 8,
      yes = str_extract(string = reference, pattern = "[0-9]{8}"),
      no = NA
    ),
    reference_doi = ifelse(
      test = n >= 9,
      yes = str_extract(string = reference, pattern = "^10.*"),
      no = NA
    ),
    reference_split = ifelse(
      test = n >= 20,
      yes = str_extract(string = reference_split, pattern = "^[A-Z].*"),
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

data_corrected$biologicalsource <-
  gsub(
    pattern = "_",
    replacement = " ",
    x =  data_corrected$biologicalsource,
    fixed = TRUE
  )

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
  mutate(organism_clean = biologicalsource
  ) %>%
  select(
    -y,
    -nchar,
    -first,
    -last,
    -capitals_x,
    -capitals_y
  ) %>%
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
    db = "coconut",
    structure_field = c("structure_inchi", "structure_smiles", "structure_name"),
    organism_field = c("organism_clean"),
    reference_field = c(
      "reference_doi",
      "reference_pubmed",
      "reference_authors",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)
