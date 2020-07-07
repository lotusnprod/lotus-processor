# title: "PharmDBK cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(data.table) #rbindlist()
library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)
library(jsonlite)

# get paths
database <- databases$get("pharmdbk")

## files
ids <- read_delim(
  file = database$sourceFiles$tsv,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(ids)[1] <- "id"
X <- ids$id

list <- list()

for (i in 1:length(X)) {
  jsonfile <- gzfile(X[i])
  
  df <- as.data.frame(fromJSON(jsonfile, simplifyDataFrame = TRUE))
  
  df2 <- df %>%
    unnest()
  
  if (nrow(df2) != 0)
    list[[i]] <- df2
}

data_original <- rbindlist(list) %>%
  data.frame()

data_original$inchi <- y_as_na(x = data_original$inchi, y = "")
data_original$smiles <- y_as_na(x = data_original$smiles, y = "")
data_original$name <- y_as_na(x = data_original$name, y = "")
data_original$name1 <- y_as_na(x = data_original$name1, y = "")

# manipulating
data_manipulated <- data_original %>%
  mutate(inchi = ifelse(
    test = !is.na(inchi),
    yes = paste("InChI=", inchi, sep = ""),
    no = inchi
  )) %>%
  select(name,
         smiles,
         inchi,
         biologicalsource = name1,
         reference_external = source)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "pha_1",
    structure_field = c("inchi", "name", "smiles"),
    reference_field = c("reference_external")
  )

# exporting
database$writeInterim(data_standard)
