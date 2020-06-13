# title: "NPCARE cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDataExternalDbSourceNpcareOriginal,
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  mutate(reference = paste(ref, ref_link, sep = "§")) %>%
  select(
    originalid = id,
    name = compounds,
    smiles = canonical_smiles,
    biologicalpart = extract,
    pubchem = pid,
    biologicalsource = species,
    reference
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npc_1",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbNpcare,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
