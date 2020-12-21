# title: "datawarrior cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("datawarrior")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# manipulating
data_manipulated <- data_original %>%
  cSplit("reference", sep = "; ") %>%
  cSplit("remarks", sep = "; ", direction = "long") %>%
  cSplit("remarks", sep = "<NL>", direction = "wide") %>%
  select(
    name,
    smiles = Smiles,
    inchi = InChI,
    biologicalsource = remarks_1,
    reference_authors = reference_1,
    reference_title = reference_2
  ) %>%
  mutate(
    biologicalsource = gsub(
      pattern = "Origin: ",
      replacement = " ",
      x = biologicalsource,
      fixed = TRUE
    )
  ) %>%
  mutate_all(as.character) %>%
  mutate_all(trimws) %>%
  mutate(
    reference_title = ifelse(
      test = reference_title == "Giftpflanzen Pflanzengifte4. ?berarbeitete Aufl.",
      yes = "Giftpflanzen Pflanzengifte, 4. ?berarbeitete Aufl.",
      no = reference_title
    )
  ) %>%
  mutate(
    reference_isbn = ifelse(
      test = reference_title == "Giftpflanzen Pflanzengifte, 4. ?berarbeitete Aufl.",
      yes = "978-3-933203-31-1",
      no = "3-8047-1053-0"
    )
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "dat_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c("reference_authors", "reference_title", "reference_isbn")
  )

# exporting
database$writeInterim(data_standard)