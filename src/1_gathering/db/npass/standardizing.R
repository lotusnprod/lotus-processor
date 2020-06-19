# title: "NPASS cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("npass")

## files
data_original_1 <- read_delim(
  file = database$sourceFiles$tsvGeneral,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_2 <- read_delim(
  file = database$sourceFiles$tsvProperties,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_3 <- read_delim(
  file = database$sourceFiles$tsvSpeciesInfo,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

data_original_4 <- read_delim(
  file = database$sourceFiles$tsvSpeciesPair,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# joining
data_original <- left_join(data_original_1, data_original_2)

data_original <- left_join(data_original, data_original_4)

data_original <- left_join(data_original, data_original_3)

# selecting
data_selected <- data_original %>%
  mutate(reference = paste(ref_id, ref_id_type, sep = "§")) %>%
  select(
    pubchem = pubchem_cid,
    np_id,
    name = pref_name,
    inchi = standard_inchi,
    standard_inchi_key,
    smiles = canonical_smiles,
    biologicalsource = org_name,
    reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_1",
    structure_field = c("name", "inchi", "smiles")
  )

data_standard$name <- y_as_na(data_standard$name, "n.a.")

# exporting
database$writeInterim(data_standard)
