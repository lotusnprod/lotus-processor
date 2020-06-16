# title: "AFROTRYP cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# files
data_original <- read_delim(
  file = pathDataExternalDbSourceAfrotrypOriginal,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  mutate(reference = paste(Reference, Publisher, sep = " ")) %>%
  select(
    uniqueid = `Compound code`,
    name = `Compound name`,
    biologicalsource = `Species name`,
    biologicalpart = `plant part`,
    reference
  )

data_selected$reference <- gsub("NA", "", data_selected$reference)

# standardizing
data_standard <-
  standardizing_original(data_selected = data_selected,
                         db = "afr_1",
                         structure_field = "name")

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbAfrotryp,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
