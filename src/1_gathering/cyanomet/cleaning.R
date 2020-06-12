# title: "CyanoMetDB cleaneR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## files
data_original <- read_delim(
  file = pathDataExternalDbSourceCyanometdb,
  delim = ",",
  escape_double = TRUE,
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

data_manipulated <- data_original %>%
  mutate(
    name = `Compound name`,
    biologicalsource = paste(ifelse(is.na(Genus),
                                    "",
                                    Genus),
                             ifelse(is.na(Species),
                                    "",
                                    Species),
                             sep = " "),
    inchi = Inchl,
    smiles = `SMILES (canonical or isomeric)`
  )

data_manipulated$reference <-
  apply(data_manipulated[, c(11, 17, 23)] , 1 , paste , collapse = "|")

data_manipulated$biologicalsource <-
  y_as_na(data_manipulated$biologicalsource, " ")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "cya_1",
    structure_field = c("name", "inchi", "smiles")
  )

# exporting
write.table(
  x = data_standard,
  file = gzfile(
    description = pathDataInterimDbCyanometdb,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
