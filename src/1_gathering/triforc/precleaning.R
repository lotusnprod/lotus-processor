#title: "TRIFORC cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "TRIFORC"
originalfile <- "0_initial_files/TriforC_original.tsv"

##paths
outpath <- paste("0_initial_files/",
                 db,
                 "_to_get.tsv",
                 sep = "")

##files
data_original <- read_delim(
  file = originalfile,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#cleaning
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

#filtering lost values to retrieve manually
data_to_get <- data_selected_final %>%
  filter(grepl("others", biologicalsource))

#exporting
write.table(
  x = data_to_get,
  file = outpath,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
