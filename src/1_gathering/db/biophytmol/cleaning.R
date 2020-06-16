# title: "Biophytmol cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceBiophytmolOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
data_selected <- data_original %>%
  select(uniqueid,
         name,
         smiles,
         biologicalsource,
         reference) %>%
  cSplit("biologicalsource", "     ") %>%
  select(uniqueid,
         name,
         smiles,
         biologicalsource = biologicalsource_1,
         reference) %>%
  mutate_all(as.character) %>%
  tibble()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "bio_2",
    structure_field = c("name", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbBiophytmol,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
