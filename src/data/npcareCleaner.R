#title: "NPCARE cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "NPCARE"
originalfile <- "0_initial_files/npcare.zip"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = originalfile,
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  mutate(reference = paste(ref, ref_link, sep = "§")) %>%
  select(
    originalid = id,
    name = compounds,
    smiles = canonical_smiles,
    biologicalpart = extract,
    pubchem = pid,
    biologicalsource = species,
    reference
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npc_1",
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
