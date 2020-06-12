#title: "MITISHAMBA cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

##files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceMitishambaOriginal),
  delim = "\t",
  trim_ws = TRUE,
  escape_backslash = TRUE
) %>%
  mutate_all(as.character)

#selecting
data_selected <- data_original %>%
  select(
    smiles,
    biologicalsource = plant_species,
    name = common_name,
    reference = authors
  )

#standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "mit_1",
    structure_field = c("name", "smiles")
  )

#exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbMitishamba,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
