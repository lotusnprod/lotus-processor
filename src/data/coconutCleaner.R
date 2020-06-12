# title: "COCONUT cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathCoconutOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
data_selected <- data_original %>%
  select(
    name,
    inchi = inchi,
    smiles = SMILES,
    biologicalsource = textTaxa,
    reference = citationDOI,
    externalDB = found_in_databases
  ) %>%
  mutate(
    biologicalsource = gsub("\\[", "", biologicalsource),
    biologicalsource = gsub("\\]", "", biologicalsource),
    biologicalsource = gsub("notax", NA, biologicalsource),
    biologicalsource = gsub("\"", "", biologicalsource),
    reference = gsub("\\[", "", reference),
    reference = gsub("\\]", "", reference),
    reference = gsub(",", "|", reference),
    reference = gsub("\"", "", reference),
    externalDB = gsub("\\[", "", externalDB),
    externalDB = gsub("\\]", "", externalDB),
    externalDB = gsub(",", "|", externalDB),
    externalDB = gsub("\"", "", externalDB),
  )

data_selected$name <- y_as_na(data_selected$name, "")
data_selected$inchi <- y_as_na(data_selected$inchi, "")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "")
data_selected$reference <- y_as_na(data_selected$reference, "")
data_selected$externalDB <- y_as_na(data_selected$externalDB, "")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "coc_1",
    structure_field = c("inchi", "smiles", "name")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathCoconutStandard,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
