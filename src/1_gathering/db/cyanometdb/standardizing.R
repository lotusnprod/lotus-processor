# title: "CyanoMetDB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("cyanometdb")

## files
data_original <- vroom(
  file = database$sourceFiles$tsv,
  delim = ",",
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

data_manipulated <- data_original %>%
  mutate(
    name = `Compound name`,
    biologicalsource = paste(ifelse(is.na(Genus),
      "",
      Genus
    ),
    ifelse(is.na(Species),
      "",
      Species
    ),
    sep = " "
    ),
    inchi = InChI,
    smiles = `SMILES (canonical or isomeric)`
  )

data_manipulated$biologicalsource <-
  y_as_na(data_manipulated$biologicalsource, " ")

data_selected <- data_manipulated %>%
  select(
    structure_name = name,
    organism_clean = biologicalsource,
    structure_inchi = inchi,
    structure_smiles = smiles,
    reference_doi_1 = `DOI_No1`,
    reference_doi_2 = `DOI_No2`,
    reference_doi_3 = `DOI_No3`,
  ) %>%
  pivot_longer(5:7) %>%
  filter(!is.na(value)) %>%
  select(structure_name,
    organism_clean,
    structure_inchi,
    structure_smiles,
    reference_doi = value
  ) %>%
  mutate(organism_clean = gsub(
    pattern = "n.a.",
    replacement = "",
    x = organism_clean,
    fixed = TRUE
  )) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "cyanometdb",
    structure_field = c("structure_name", "structure_inchi", "structure_smiles"),
    organism_field = "organism_clean",
    reference_field = c("reference_doi")
  )

# exporting
database$writeInterim(data_standard)
