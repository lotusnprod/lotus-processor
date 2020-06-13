#title: "STREPTOMEDB cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceStreptomedbCompiled),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(uniqueid,
         name,
         #inchi = ,
         #inchikey = ,
         smiles,
         #cas = ,
         pubchem,
         reference = pubmedid,
         biologicalsource)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "str_1",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbStreptomedb,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
