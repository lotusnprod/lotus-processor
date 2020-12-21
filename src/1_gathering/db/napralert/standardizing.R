# title: "NAPRALERT cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(Hmisc)

# get paths
database <- databases$get("napralert")

## files
dataOriginal <- read_delim(
  file = gzfile(database$sourceFiles$tsvOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    biologicalsource = paste(
      capitalize(string = tolower(`?genus`)),
      tolower(`?species`),
      capitalize(string = tolower(`?family`))
    ),
    inchi = NA,
    reference_doi = NA
  ) %>%
  select(
    name = `?compound_name`,
    inchi,
    biologicalsource,
    reference_title = `?title`,
    reference_authors = `?authors`,
    reference_doi,
    reference_journal = `?journal`
  ) %>%
  mutate_all(as.character)

dataMatched <- read_delim(
  file = gzfile(database$sourceFiles$tsvMatched),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(
    name = NA,
    reference_title = NA,
    reference_authors = NA,
    reference_journal = NA,
  ) %>%
  select(
    name,
    inchi = InChI,
    biologicalsource = TaxonName,
    reference_title,
    reference_authors,
    reference_journal,
    reference_doi = DOI
  ) %>%
  mutate_all(as.character)

# manipulating
dataJoined <- bind_rows(dataMatched, dataOriginal)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = dataJoined,
    db = "nap_1",
    structure_field = c("name", "inchi"),
    reference_field = c(
      "reference_doi",
      "reference_authors",
      "reference_title",
      "reference_journal"
    )
  )

# exporting
database$writeInterim(data_standard)