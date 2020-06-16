# title:"NPEDIA cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceNpediaOriginal),
  delim = "\t",
  escape_double = TRUE,
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(
    uniqueid = ID,
    name = Name,
    biologicalsource = Source,
    inchi = InChI,
    cas = `CAS No.`,
    reference = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npe_1",
    structure_field = c("name", "inchi")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbNpedia,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
