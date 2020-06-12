#title: "PLANTCYC cleaneR"

#loading
##functions
source("../../functions.R")


##db
db <- "PLANTCYC"

##paths
filenames <- list.files("0_initial_files/",
                        pattern = "*.tsv",
                        full.names = TRUE)

outpath <- "PLANTCYC_std.tsv.zip"

data_standard <- do.call("rbind",
                         lapply(filenames,
                                function(x) {
                                  dat <- read.csv(x,
                                                  header = TRUE,
                                                  sep = "\t")
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
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
