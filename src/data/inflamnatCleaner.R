#title: "INFLAMNAT cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "INFLAMNAT"
originalfile <- "0_initial_files/ci8b00560_si_001.xlsx"

##paths
inpath <- paste("0_initial_files/",
                originalfile,
                sep = "")

outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_excel(originalfile,
                            sheet = 1) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  select(
    uniqueid = Index,
    name = Name,
    smiles = SMILES,
    biologicalsource = Origin,
    pubchem = CID,
    reference = Reference
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "inf_1",
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
