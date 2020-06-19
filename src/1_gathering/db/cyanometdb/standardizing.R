# title: "CyanoMetDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("cyanometdb")

## files
data_original <- read_delim(
  file = database$sourceFiles$tsv,
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
database$writeInterim(data_standard)
