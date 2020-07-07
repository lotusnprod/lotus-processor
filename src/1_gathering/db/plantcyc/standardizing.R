# title: "PLANTCYC cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("plantcyc")

data_standard <- do.call("rbind",
                         lapply(pathDataExternalDbSourcePlantcycOriginal,
                                function(x) {
                                  dat <- read_delim(
                                    file = gzfile(x),
                                    delim = "\t",
                                    escape_double = FALSE,
                                    trim_ws = TRUE
                                  )
                                  dat
                                }))

data_standard$name <- NA
data_standard$reference_external <- "PLANTCYC"

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_standard,
    db = "pla_1",
    structure_field = c("name", "inchi"),
    reference_field = c("reference_external")
  )

#exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbPlantcyc,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
