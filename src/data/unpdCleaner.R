#title: "UNPD cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "UNPD"

##paths
originalfile <- "0_initial_files/unpd_translated.tsv.zip"

outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = gzfile(originalfile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

#selecting
##atomizing references
data_selected <- data_original %>%
  mutate(reference = gsub("(\\(\\d+).\\s","| ",ref),
         reference = sub("\\| ","",reference)) %>% 
  select(
    biologicalsource = ln_reduced,
    reference,
    inchi = InChI,
    smiles = SMILES
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "unp_1",
    structure_field = c("inchi", "name", "smiles")
  )

#exporting
write.table(
  x = data_standard,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
