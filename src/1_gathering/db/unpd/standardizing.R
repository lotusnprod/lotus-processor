# title: "UNPD cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("unpd")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceUnpdIntegrated),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
## atomizing references
data_selected <- data_original %>%
  mutate(reference = gsub("(\\(\\d+).\\s", "|", ref),) %>%
  cSplit("reference", sep = "|", direction = "long") %>%
  mutate_all(as.character) %>%
  select(
    biologicalsource = ln_reduced,
    reference_unsplittable = reference,
    inchi = InChI,
    smiles = SMILES
  ) %>%
  data.frame()

data_manipulated <- data_selected %>%
  mutate(
    reference_external = ifelse(
      test = reference_unsplittable == "Retrieved from CNPD",
      yes = reference_unsplittable,
      no = NA
    ),
    reference_unsplittable = gsub(
      pattern = "Retrieved from CNPD",
      replacement = "",
      x = reference_unsplittable,
      fixed = TRUE
    )
  )

data_manipulated$reference_unsplittable <-
  y_as_na(x = data_manipulated$reference_unsplittable, y = "")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "unp_1",
    structure_field = c("inchi", "name", "smiles"),
    reference_field = c("reference_unsplittable", "reference_external")
  )

# exporting
database$writeInterim(data_standard)
