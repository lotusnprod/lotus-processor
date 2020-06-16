# title: "UNPD cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDataExternalDbSourceUnpdOriginal_1,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)

data_original_ISDB <- read_delim(
  file = pathDataExternalDbSourceUnpdOriginal_2,
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
    description = pathDataExternalDbSourceUnpdCompiled,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
