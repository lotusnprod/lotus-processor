#title: "NANPDB cleaneR"

#loading
##functions
source("../../functions.R")

##db
db <- "NANPDB"
originalfile <- "0_initial_files/NANPDB_scraped.tsv.zip"

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
  mutate_all(as.character) %>%
  mutate(id = row.names(.))

#selecting
data_selected <- data_original %>%
  mutate(name = NA) %>%
  select(
    uniqueid = id,
    name,
    smiles = SMILES.,
    biologicalsource = Source.Species.Information,
    pubchem = PubChem.,
    reference = Reference.information
  ) %>%
  cSplit("biologicalsource",
         " Known use:",
         stripWhite = FALSE,
         fixed = FALSE) %>%
  select(uniqueid,
         name,
         smiles,
         biologicalsource = biologicalsource_1,
         pubchem,
         reference) %>%
  mutate(biologicalsource = gsub("Source: ", "", biologicalsource)) %>%
  mutate_all(as.character) %>%
  cSplit("reference",
         "Title: ",
         stripWhite = FALSE,
         fixed = FALSE) %>%
  cSplit("reference_2",
         " PubMed: ",
         stripWhite = FALSE,
         fixed = FALSE) %>%
  cSplit("reference_1",
         " Reference: ",
         stripWhite = FALSE,
         fixed = FALSE) %>%
  mutate(reference = paste(reference_2_1, reference_1_2, sep = ", ")) %>%
  tibble()

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "nan_1",
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
