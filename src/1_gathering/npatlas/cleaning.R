# title: "NPATLAS cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDataExternalDbSourceNpatlasOriginal,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# selecting
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

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "npa_2",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(description = pathDataInterimDbNpatlas,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
