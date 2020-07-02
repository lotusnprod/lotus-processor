# title: "TipDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("carotenoiddb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# manipulating
data_manipulated <- data_original %>%
  select(
    uniqueid = Entry,
    name = `Carotenoid Name`,
    inchi = InChI,
    inchikey = InChIKey,
    smiles = `Canonical SMILES`,
    cas = CAS,
    biologicalsource = `Source organisms`,
    reference = References
  )

data_manipulated$biologicalsource <-
  gsub("(\\))([A-Z])",
       ")SPLIT\\2",
       data_manipulated$biologicalsource)

data_manipulated <- data_manipulated %>%
  cSplit(
    "biologicalsource",
    sep = "SPLIT",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  mutate(across(.cols = everything(), as.character))

data_manipulated_long <- data_manipulated %>%
  gather(c(8:ncol(.)),
         key = "n",
         value = "biologicalsource") %>%
  group_by(uniqueid) %>%
  select(-n) %>%
  distinct(biologicalsource, .keep_all = TRUE) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(biologicalsource) | !n > 1) %>%
  select(-n) %>%
  arrange(uniqueid)

data_manipulated_long$reference <-
  gsub("(Ref.\\d* : )", "SPLIT\\1", data_manipulated_long$reference)

data_manipulated_long_ref <- data_manipulated_long %>%
  cSplit("reference",
         sep = "SPLIT",
         stripWhite = FALSE,
         fixed = FALSE) %>%
  mutate_all(as.character) %>%
  gather(c(8:ncol(.)),
         key = "n",
         value = "reference") %>%
  select(-n) %>%
  group_by(uniqueid) %>%
  distinct(biologicalsource, reference, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(uniqueid, biologicalsource) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(reference) | !n > 1) %>%
  select(-n) %>%
  arrange(uniqueid)

data_manipulated_long_ref_unique <- data_manipulated_long_ref %>%
  mutate(
    refnum = str_extract(data_manipulated_long_ref$reference, "(Ref.\\d*)"),
    biorefnum = str_extract(data_manipulated_long_ref$biologicalsource, "(Ref.\\d*)")
  ) %>%
  filter(refnum == biorefnum | is.na(refnum))

data_manipulated_long_ref_unique$reference <-
  gsub("(Ref.\\d* : )",
       "",
       data_manipulated_long_ref_unique$reference)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated_long_ref_unique,
    db = "car_1",
    structure_field = c("inchi", "name", "smiles")
  )

# exporting
database$writeInterim(data_standard)
