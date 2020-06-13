# title: "PLANTCYC cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

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

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_standard,
    db = "pla_1",
    structure_field = c("name", "inchi")
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
