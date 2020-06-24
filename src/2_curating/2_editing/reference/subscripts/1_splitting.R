# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions/helpers.R")
source("functions/reference.R")

library(dplyr)
library(readr)
library(stringr)
library(splitstackshape)
library(tidyr)

# loading files
## reference
dataReference <- read_delim(
  file = gzfile(pathOriginalRef),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# splitting when MULTIPLE REFERENCES (|)
dataReferenceSplit <- dataReference %>%
  mutate(newreference = referenceOriginal) %>%
  cSplit(splitCols = "newreference",
         sep = "|") %>%
  mutate_all(as.character) %>%
  tibble()

# pivoting
dataReferenceLong <- dataReferenceSplit %>%
  rowwise() %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  filter(!is.na(value)) %>%
  ungroup()

dataReferenceLong$value <- y_as_na(x = dataReferenceLong$value,
                                   y = "")

dataReferenceLong$value <- y_as_na(x = dataReferenceLong$value,
                                   y = "NA")

dataReferenceLong <- dataReferenceLong %>%
  filter(!is.na(value))

# splitting when multiple fields FOR A SINGLE REFERENCE (§)
dataReferenceLongSplit <- dataReferenceLong %>%
  cSplit(splitCols = "value",
         sep = "§") %>%
  rowwise() %>%
  pivot_longer(
    cols = (ncol(.) - 1):ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  ungroup()

# selecting fields probably corresponding to a title or authors list (human)
dataReferenceFillAuto <- dataReferenceLongSplit %>%
  filter(
    !grepl(pattern = "^KNApSAcK Database$",
           x = value) &
      !grepl(pattern = "^DrDuke$",
             x = value) &
      !grepl(pattern = "^patent:",
             x = value) &
      !grepl(pattern = "^CHEBI:",
             x = value) &
      !grepl(pattern = "^doi:",
             x = value) &
      !grepl(pattern = "^10\\.\\d{4,9}",
             x = value)  &
      !grepl(pattern = "^pubmed:",
             x = value) &
      !grepl(pattern = "^PUBMED ",
             x = value) &
      !grepl(pattern = "^[0-9]*$",
             x = value) &
      !grepl(pattern = "^\\d+(;\\d+)*$",
             x = value)  &
      !grepl(pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
             x = value)
  ) %>% 
  mutate_all(as.character)

# selecting fields probably corresponding to pubmed ID
dataReferenceFillPubmed <- dataReferenceLongSplit %>%
  filter(
    grepl(pattern = "^pubmed:",
          x = value) |
      grepl(pattern = "^PUBMED ",
            x = value) |
      grepl(pattern = "^[0-9]*$",
            x = value) |
      grepl(pattern = "^\\d+(;\\d+)*$",
            x = value) |
      grepl(pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
            x = value)
  ) %>%
  cSplit(splitCols = "value",
         sep = ";") %>%
  select(-level) %>%
  mutate_all(as.character) %>%
  group_by(referenceOriginal) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "taxonomy",
    values_drop_na = TRUE
  ) %>%
  ungroup() %>%
  mutate(value = gsub(
    pattern = "pubmed:",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "PUBMED ",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "http://www.ncbi.nlm.nih.gov/pubmed/",
    replacement = "",
    x = value
  )) %>%
  mutate(value = gsub(
    pattern = "?term=",
    replacement = "",
    x = value
  )) %>% 
  mutate_all(as.character)

# selecting fields corresponding to DOIs
dataReferenceFillDoi <- dataReferenceLongSplit %>%
  filter(grepl("^doi:", value) |
           grepl("^10\\.\\d{4,9}", value)) %>% 
  mutate_all(as.character)
