# title: "CyanoMetDB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(data.table)
library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("cyanometdb")

## files
data_original <- data.table::fread(
  file = database$sourceFiles$tsv,
  strip.white = FALSE,
  encoding = "Latin-1"
) %>%
  dplyr::mutate_all(as.character)

data_manipulated <- data_original |>
  dplyr::mutate(
    name = IUPAC_name,
    biologicalsource = paste(
      ifelse(is.na(Genus),
        "",
        Genus
      ),
      ifelse(is.na(Species),
        "",
        Species
      ),
      sep = " "
    )
  )

data_manipulated$biologicalsource <-
  y_as_na(data_manipulated$biologicalsource, " ")

data_selected <- data_manipulated |>
  dplyr::select(
    structure_name = name,
    organism_clean = biologicalsource,
    structure_inchi = InChI,
    structure_smiles = SMILES,
    reference_doi_1 = `DOI_No1`,
    reference_doi_2 = `DOI_No2`,
    reference_doi_3 = `DOI_No3`,
  ) |>
  tidyr::pivot_longer(5:7) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::select(structure_name,
    organism_clean,
    structure_inchi,
    structure_smiles,
    reference_doi = value
  ) |>
  dplyr::mutate(organism_clean = gsub(
    pattern = "n.a.",
    replacement = "",
    x = organism_clean,
    fixed = TRUE
  )) |>
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "cyanometdb",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c("reference_doi")
  )

# exporting
database$writeInterim(data_standard)
