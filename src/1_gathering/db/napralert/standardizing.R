# title: "NAPRALERT cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("napralert")

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
    inchi = InChI,
    inchikey = InChIKey,
    biologicalsource = TaxonName,
    reference_doi = DOI
  ) %>%
  data.frame()

data_manipulated$name <- NA

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "nap_1",
    structure_field = c("inchi"),
    reference_field = c("reference_doi")
  )

# exporting
database$writeInterim(data_standard)
