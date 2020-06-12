# title: "DNP cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDnpOriginal,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

## selecting
data_selected <- data_original %>%
  select(
    uniqueid = CRC_Number,
    name = Molecule_Name,
    inchi = MolfileName,
    biologicalsource = Biological_Source
  ) %>%
  distinct(uniqueid, .keep_all = TRUE)

## standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "dnp_1",
    structure_field = c("name", "inchi")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDnpStandard,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
