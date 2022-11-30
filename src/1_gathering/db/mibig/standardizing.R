# title: "MIBIG cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(jsonlite)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("mibig")

fileInZip <-
  function(inZip) {
    outFile <- list()
    fileList <- unzip(inZip, list = TRUE)
    fileList <- fileList |>
      filter(Length != 0)
    for (i in seq_along(fileList$Name)) {
      if (grepl(".json", fileList[i, 1])) {
        oFa <- fromJSON(unz(inZip, fileList[i, 1]))
        outFile[[i]] <- oFa
      }
    }
    return(outFile)
  }

df <- fileInZip(inZip = database$sourceFiles$data)

x <- seq_along(df)

getid <- function(x) {
  j <- df[[x]][["cluster"]][["mibig_accession"]]
  j
}

id <- lapply(
  FUN = getid,
  X = x
)

getsmiles <- function(x) {
  j <- df[[x]][["cluster"]][["compounds"]][["chem_struct"]]
  j
}

smiles <- lapply(
  FUN = getsmiles,
  X = x
)

getorganism <- function(x) {
  j <- df[[x]][["cluster"]][["organism_name"]]
  j
}

organism <- lapply(
  FUN = getorganism,
  X = x
)

getreference <- function(x) {
  j <- df[[x]][["cluster"]][["publications"]]
  j
}

reference <- lapply(
  FUN = getreference,
  X = x
)

data <- tibble(id, smiles, organism, reference) %>%
  filter(smiles != "NULL") %>%
  unnest(id) %>%
  unnest(smiles) %>%
  unnest(organism) %>%
  unnest(reference) %>%
  distinct(id, .keep_all = TRUE) %>%
  filter(!is.na(smiles)) %>%
  select(id,
    structure_smiles = smiles,
    organism_clean = organism,
    reference
  ) %>%
  mutate(
    reference_doi = gsub(
      pattern = "doi:",
      replacement = "",
      x = str_extract(string = reference, pattern = "^doi.*")
    ),
    reference_pubmed = gsub(
      pattern = "pubmed:",
      replacement = "",
      x = str_extract(string = reference, pattern = "^pubmed.*")
    ),
    reference_external = ifelse(
      test = is.na(reference_doi) & is.na(reference_pubmed),
      yes = reference,
      no = NA
    )
  )

data$name <- NA

# standardizing
data_standard <- standardizing_original(
  data_selected = data,
  db = "mibig",
  structure_field = "structure_smiles",
  organism_field = "organism_clean",
  reference_field = c("reference_doi", "reference_pubmed", "reference_external")
)

# exporting
database$writeInterim(data_standard)
