#title: "PAMDB cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "PAMDB"

##paths
inpath <- "0_initial_files/PAMDB_scraped.tsv.zip"
  
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = gzfile(inpath),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

#replacing Not available
data_original[] <-
  lapply(data_original, function(x)
    gsub("Not Available", NA, x))

#selecting
data_selected <- data_original %>%
  select(
    uniqueid = MetaboliteID,
    name = Name,
    inchi = InChI,
    smiles = SMILES,
    cas = CAS_number
  ) %>%
  mutate(biologicalsource = "Pseudomonas aeruginosa")

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "pam_1",
    structure_field = c("inchi", "smiles", "cas")
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
