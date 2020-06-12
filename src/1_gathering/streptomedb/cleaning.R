#title: "STREPTOMEDB cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "STREPTOMEDB"
originalfile <- "0_initial_files/streptomedb.tsv.zip"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = gzfile(originalfile),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#applying
##selecting
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

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "str_1",
    structure_field = c("name", "smiles")
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
