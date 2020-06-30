# title: "METABOLIGHTS studies cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(jsonlite)
library(pbmcapply)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("metabolights")

files <-
  dir(path = pathDataExternalDbSourceMetabolightsStudiesScrapedDir,
      pattern = "*.json")

filenames <-
  list.files(path = pathDataExternalDbSourceMetabolightsStudiesScrapedDir,
             pattern = "*.json",
             full.names = TRUE)

# just to get problematic entries

# for (i in files){
# i <- stream_in(file(i), verbose = FALSE)}

df <- rbind(lapply(filenames, fromJSON))

x <- 1:length(df)

getid <- function(x) {
  j <- df[[x]][["id"]]
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

getname <- function(x) {
  j <- df[[x]][["name"]]
  j
}

name <- pbmclapply(
  FUN = getname,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
)

getinchi <- function(x) {
  j <- df[[x]][["inchi"]]
  j
}

inchi <- pbmclapply(
  FUN = getinchi,
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
  j <- df[[x]][["smiles"]]
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

getspecies <- function(x) {
  j <- names(df[[1]][["species"]])
  j
}

species <- pbmclapply(
  FUN = getspecies,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE, 
  ignore.interactive = TRUE
)

data <- tibble(id, name, inchi, smiles, species) %>%
  unnest(species) %>%
  mutate_all(as.character)

species <- data %>% distinct(species)

write.table(
  x = data,
  file = gzfile(
    description = database$sourceFiles$tsvStudies,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
