# title: "NAPRALERT cleaneR"

# loading paths
source("paths.R")
source("r/capitalize.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("napralert")

## files
dataOriginal <- readr::read_delim(
  file = gzfile(description = database$sourceFiles$tsvOriginal),
  col_types = cols(.default = "c"),
  escape_backslash = TRUE,
  delim = "\t"
) |>
  dplyr::mutate(
    biologicalsource = paste(
      capitalize(string = tolower(`?genus`)),
      tolower(`?species`),
      capitalize(string = tolower(`?family`))
    ),
    inchi = NA,
    reference_doi = NA
  ) |>
  dplyr::select(
    structure_name = `?compound_name`,
    structure_inchi = inchi,
    organism_clean = biologicalsource,
    reference_title = `?title`,
    reference_authors = `?authors`,
    reference_doi,
    reference_journal = `?journal`
  )

dataMatched <- readr::read_delim(
  file = gzfile(description = database$sourceFiles$tsvMatched),
  col_types = cols(.default = "c")
) |>
  dplyr::mutate(
    name = NA,
    reference_title = NA,
    reference_authors = NA,
    reference_journal = NA,
  ) |>
  dplyr::select(
    structure_name = name,
    structure_inchi = InChI,
    organism_clean = TaxonName,
    reference_title,
    reference_authors,
    reference_journal,
    reference_doi = DOI
  )

# manipulating
dataJoined <- dplyr::bind_rows(dataMatched, dataOriginal)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = dataJoined,
    db = "napralert",
    structure_field = c("structure_name", "structure_inchi"),
    organism_field = "organism_clean",
    reference_field = c(
      "reference_doi",
      "reference_authors",
      "reference_title",
      "reference_journal"
    )
  )

# exporting
database$writeInterim(data_standard)
