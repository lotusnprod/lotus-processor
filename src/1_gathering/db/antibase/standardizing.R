# title: "AntiBase cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)

# get paths
database <- databases$get("antibase")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv)
)

data_smiles <- read_delim(
  file = gzfile(database$sourceFiles$smi)
)

# selecting
data_selected <- data_original %>%
  select(
    MOL_ID,
    name = NAME,
    biologicalsource = SOURCE,
    reference_publishingDetails = REFERENCES
  ) %>%
  mutate(biologicalsource = gsub(
    pattern = "\\[+[[:alpha:]]+\\]",
    replacement = "",
    x = biologicalsource
  )) %>%
  cSplit(
    "biologicalsource",
    sep = ";",
    direction = "long",
    fixed = TRUE
  ) %>%
  cSplit(
    "biologicalsource",
    sep = ",",
    direction = "long",
    fixed = TRUE
  )

data_corrected <- data_selected %>%
  left_join(data_smiles) %>%
  select(structure_name = name,
         organism_clean = biologicalsource,
         organism_dirty = biologicalsource,
         reference_publishingDetails,
         structure_smiles = SMILES) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected,
    db = "antibase",
    structure_field = c("structure_smiles", "structure_name"),
    organism_field = c("organism_clean", "organism_dirty"),
    reference_field = c(
      "reference_publishingDetails"
    )
  )

# exporting
database$writeInterim(data_standard)
