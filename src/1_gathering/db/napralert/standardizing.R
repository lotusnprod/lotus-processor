# title: "NAPRALERT cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(Hmisc, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

# get paths
database <- databases$get("napralert")

## files
dataOriginal <- read_delim(
  file = gzfile(database$sourceFiles$tsvOriginal),
  delim = "\t",
  col_types = cols(.default = "c")
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
  )

dataMatched <- vroom(
  file = gzfile(database$sourceFiles$tsvMatched),
  delim = "\t",
  col_types = cols(.default = "c")
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
  )

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
