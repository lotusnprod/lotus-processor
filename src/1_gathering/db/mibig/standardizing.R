# title: "MIBIG cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")
source("functions/parallel.R")

library(dplyr)
library(jsonlite)
library(pbmcapply)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("mibig")

df <- rbind(lapply(database$sourceFiles$tsv, fromJSON))

x <- 1:length(df)

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
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
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
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
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
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
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
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
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
         smiles,
         biologicalsource = organism,
         reference)

data$name <- NA

# standardizing
data_standard <- standardizing_original(
  data_selected = data,
  db = "mib_1",
  structure_field = c("name", "smiles")
)

# exporting
database$writeInterim(data_standard)
