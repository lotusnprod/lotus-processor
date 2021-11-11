# title: "MIBIG cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")
source("r/parallel.R")

library(dplyr)
library(jsonlite)
library(parallel)
library(pbmcapply)
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

id <- pbmclapply(
  FUN = getid,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = numCores,
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE,
  ignore.interactive = TRUE,
  mc.style = "txt",
  mc.substyle = 1
)

getsmiles <- function(x) {
  j <- df[[x]][["cluster"]][["compounds"]][["chem_struct"]]
  j
}

smiles <- pbmclapply(
  FUN = getsmiles,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = numCores,
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE,
  ignore.interactive = TRUE,
  mc.style = "txt",
  mc.substyle = 1
)

getorganism <- function(x) {
  j <- df[[x]][["cluster"]][["organism_name"]]
  j
}

organism <- pbmclapply(
  FUN = getorganism,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = numCores,
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE,
  ignore.interactive = TRUE,
  mc.style = "txt",
  mc.substyle = 1
)

getreference <- function(x) {
  j <- df[[x]][["cluster"]][["publications"]]
  j
}

reference <- pbmclapply(
  FUN = getreference,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = numCores,
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE,
  ignore.interactive = TRUE,
  mc.style = "txt",
  mc.substyle = 1
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
