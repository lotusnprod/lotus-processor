# title: "TRIFORC precleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDataExternalDbSourceTriforcOriginal,
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
