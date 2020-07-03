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

for (i in 1:length(X)){
  
jsonfile <- gzfile(X[i])

df <- as.data.frame(fromJSON(jsonfile, simplifyDataFrame = TRUE))

df2 <- df %>% 
  unnest()

if (nrow(df2)!=0)
list[[i]] <- df2
}

data_original <- rbindlist(list) %>% 
  data.frame()

# manipulating
data_manipulated <- data_original %>%
  select(
    name,
    smiles,
    inchi,
    biologicalsource = name1,
    reference = source
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "pha_1",
    structure_field = c("inchi", "name", "smiles")
  )

# exporting
database$writeInterim(data_standard)
