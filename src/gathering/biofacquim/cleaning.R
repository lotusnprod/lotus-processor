# title: "BIOFACQUIM cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathBiofacquimOriginal,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    name = Name,
    smiles = SMILES,
    biologicalsource = Specie,
    reference = Reference
  )

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_1",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathBiofacquimStandard,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
