# title: "UNPD cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("unpd")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsvJo,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

data_original_ISDB <- read_delim(
  file = database$sourceFiles$tsvPm,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
data_translated <-
  left_join(data_original,
            data_original_ISDB,
            by = c("inchik" = "inchik")) %>%
  select(unpdid,
         inchik,
         ln_reduced,
         ref,
         InChI,
         SMILES)

data_translated[] <-
  lapply(data_translated, function(x)
    gsub("\r\n", " ", x))

data_translated[] <-
  lapply(data_translated, function(x)
    gsub("\r", " ", x))

data_translated[] <-
  lapply(data_translated, function(x)
    gsub("\n", " ", x))

# exporting
write.table(
  x = data_translated,
  file = gzfile(
    description = pathDataExternalDbSourceUnpdIntegrated,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)