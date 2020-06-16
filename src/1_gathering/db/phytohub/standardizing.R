# title: "Phytohub cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourcePhytohubOriginal),
  delim = "\t",
  escape_double = TRUE,
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_original,
    db = "phy_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbPhytohub,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
