#title: "UNPD cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "UNPD"

##paths
originalfile <- "0_initial_files/unpd_final.csv.zip"

originalfileISDB <- "0_initial_files/UNPD_DB.csv.zip"

outpath <- "0_initial_files/unpd_translated.tsv.zip"

##files
data_original <- read_delim(
  file = originalfile,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

data_original_ISDB <- read_delim(
  file = originalfileISDB,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

#selecting
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

#exporting
write.table(
  x = data_translated,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
