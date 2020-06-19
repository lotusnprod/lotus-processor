# title: "TRIFORC precleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("triforc")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv1,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# cleaning
data_selected_final <-  data_original %>%
  select(
    name = Name,
    cas = CAS,
    pubchem = `PubChem CID`,
    biologicalsource = Plant
  ) %>%
  cSplit("biologicalsource", ", ") %>%
  pivot_longer(
    4:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "biologicalsource",
    values_drop_na = TRUE
  )

# filtering lost values to retrieve manually
dataToGet <- data_selected_final %>%
  filter(grepl("others", biologicalsource))

# exporting
write.table(
  x = dataToGet,
  file = pathDataExternalDbSourceTriforcToGet,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
