# title: "UNPD cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)

# get paths
database <- databases$get("unpd")

## files
data_original <-
  readr::read_delim(file = database$sourceFiles$tsvJo)

data_original_ISDB <-
  readr::read_delim(file = database$sourceFiles$tsvPm)

# selecting
data_translated <-
  dplyr::left_join(
    data_original,
    data_original_ISDB,
    by = c("inchik" = "inchik")
  ) |>
  dplyr::select(
    unpdid,
    inchik,
    ln_reduced,
    ref,
    InChI,
    SMILES
  )

data_translated[] <-
  lapply(data_translated, function(x) {
    gsub(
      pattern = "\r\n",
      replacement = " ",
      x = x
    )
  })

data_translated[] <-
  lapply(data_translated, function(x) {
    gsub(
      pattern = "\r",
      replacement = " ",
      x = x
    )
  })

data_translated[] <-
  lapply(data_translated, function(x) {
    gsub(
      pattern = "\n",
      replacement = " ",
      x = x
    )
  })

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
