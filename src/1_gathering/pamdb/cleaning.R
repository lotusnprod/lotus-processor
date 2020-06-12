#title: "PAMDB cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "PAMDB"
originalfile <- "0_initial_files/PaMet.xlsx"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_excel(originalfile) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  select(
    uniqueid = MetID,
    name = Name,
    inchi = InChI,
    smiles = SMILES,
    cas = `CAS number`,
    reference = References
  ) %>%
  mutate(biologicalsource = "Pseudomonas aeruginosa")

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "pam_1",
    structure_field = c("name", "inchi", "smiles")
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
