# title: "TMMC cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(readxl)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("tmmc")

## files
data_original <- read_excel(database$sourceFiles$tsv,
                            sheet = 1) %>%
  mutate_all(as.character)

data_original_long <- data_original %>%
  cSplit("CSID", "|") %>%
  pivot_longer(17:ncol(.)) %>%
  filter(!is.na(value)) %>%
  distinct(SCIENCE, COMPOUND, .keep_all = TRUE) %>%
  mutate(
    name = COMPOUND,
    biologicalsource = str_extract(SCIENCE, "(?<=\\[).+?(?=\\])"),
    reference = LINK
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(
    biologicalsource = gsub("<i>", "", biologicalsource),
    biologicalsource = gsub("</i>", "", biologicalsource)
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original_long,
    db = "tmm_1",
    structure_field = c("name")
  )

# exporting
database$writeInterim(data_standard)
