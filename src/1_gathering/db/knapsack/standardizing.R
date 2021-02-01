# title: "KNAPSACK cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("knapsack")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t"
) %>%
  mutate_all(as.character)

## applying
data_selected <- data_original %>%
  select(
    structure_name = Name,
    uniqueid = C_ID,
    structure_inchi = InChICode,
    structure_smiles = SMILES,
    organism_clean = Organism,
    reference_original = Reference
  ) %>%
  mutate(reference_split = ifelse(
    test = grepl(
      pattern = ".*et al",
      x = reference_original
    ),
    yes =
      trimws(x = sub(
        pattern = "^ /",
        replacement = "",
        x = sub(
          pattern = "^,",
          replacement = "",
          x = sub(
            pattern = "^\\.",
            replacement = "",
            x = sub(
              pattern = ".*et al",
              replacement = "",
              x = reference_original
            )
          )
        )
      )),
    no = NA
  )) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "kna_1",
    structure_field = c("structure_name", "structure_inchi", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c("reference_original", "reference_split")
  )

# exporting
database$writeInterim(data_standard)
