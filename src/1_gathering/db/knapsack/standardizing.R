# title: "KNAPSACK cleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceKnapsackOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
) %>%
  mutate_all(as.character)

## applying
data_selected <- data_original %>%
  select(
    name = Name,
    uniqueid = C_ID,
    inchi = InChICode,
    smiles = SMILES,
    biologicalsource = Organism,
    reference = Reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "kna_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbKnapsack,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
