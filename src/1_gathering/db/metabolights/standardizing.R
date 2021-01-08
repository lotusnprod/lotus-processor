# title: "Metabolights cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

groundhog.library(Hmisc, date = groundhog.day) # neeed for capitalize()
groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

# get paths
database <- databases$get("metabolights")

data_clean_final <- vroom(
  file = gzfile(database$sourceFiles$tsvPrecleaned),
  delim = "\t"
)

species_studies <- vroom(
  file = gzfile(database$sourceFiles$tsvStudies),
  delim = "\t"
) %>%
  distinct(species)

species <- capitalize(species_studies$species)

data_selected <- data_clean_final %>%
  filter(!biologicalsource %in% species)

data_selected$inchi <- y_as_na(data_selected$inchi, "NULL")
data_selected$name <- y_as_na(data_selected$name, "NULL")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "NULL")
data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "reference compound")

data_selected <- data_selected %>%
  filter(!is.na(biologicalsource)) %>%
  mutate(reference_external = reference) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "met_1",
    structure_field = c("name", "inchi"),
    reference_field = c("reference_external")
  )

# exporting
database$writeInterim(data_standard)
