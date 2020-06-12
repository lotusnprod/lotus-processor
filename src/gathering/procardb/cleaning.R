#title: "PROCARDB cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "PROCARDB"
originalfile <- "PROCARDB_scraped.tsv.zip"

##paths
inpath <- paste("0_initial_files/",
                originalfile,
                sep = "")

outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = gzfile(inpath),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#selecting
##atomizing ref
data_selected <- data_original %>%
  mutate(
    reference = gsub("(\\d+\\.)([[:alpha:]])", "| \\2", REFERENCES),
    reference = gsub("(\\d+\\.)(\\s)([[:alpha:]])", "| \\3", reference),
    reference = gsub("(\\d+\\.)(\\s)(\\s)([[:alpha:]])", "| \\4", reference),
    reference = sub("\\| ", "", reference)
  ) %>%
  select(
    uniqueid = column_label,
    name = `CAROTENOID NAME`,
    biologicalsource,
    #inchi = InChI #is an inchikey!!!
    smiles = `CANONICAL SMILES`,
    pubchem = `PUBCHEM ID`,
    reference
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "pro_1",
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
