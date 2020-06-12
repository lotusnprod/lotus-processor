#title: "NPATLAS cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "NPATLAS"
originalfile <- "0_initial_files/np_atlas_2019_12.tsv"

##paths
outpath <- paste(db,
                 "_std.tsv.zip",
                 sep = "")

##files
data_original <- read_delim(
  file = originalfile,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  mutate(organism = paste(Genus,
                          `Origin Species`,
                          sep = " "),
         reference = `Isolation Reference DOI`) %>%
  select(
    NPAID,
    name = Names,
    InChIKey,
    inchi = InChI,
    smiles = SMILES,
    biologicalsource = organism,
    reference
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_2",
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
