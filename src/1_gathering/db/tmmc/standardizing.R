# title: "TMMC cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(readxl)
library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("tmmc")

## files
data_original <- read_excel(database$sourceFiles$tsv,
  sheet = 1
) %>%
  mutate_all(as.character)

data_original_long <- data_original %>%
  cSplit("CSID", "|") %>%
  pivot_longer(17:ncol(.)) %>%
  filter(!is.na(value)) %>%
  distinct(SCIENCE, COMPOUND, .keep_all = TRUE) %>%
  mutate(
    name = COMPOUND,
    biologicalsource = str_extract(SCIENCE, "(?<=\\[).+?(?=\\])"),
    reference_pubmed = str_extract(string = LINK, pattern = "[0-9]{6,9}")
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(
    biologicalsource = gsub("<i>", "", biologicalsource),
    biologicalsource = gsub("</i>", "", biologicalsource)
  ) %>%
  data.frame() %>%
  select(
    structure_name = name,
    organism_clean = biologicalsource,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original_long,
    db = "tmm_1",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = "reference_pubmed"
  )

# exporting
database$writeInterim(data_standard)
