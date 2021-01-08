# title: "CyanoMetDB cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

# get paths
database <- databases$get("cyanometdb")

## files
data_original <- vroom(
  file = database$sourceFiles$tsv,
  delim = ",",
  trim_ws = FALSE
) %>%
  mutate_all(as.character)

data_manipulated <- data_original %>%
  mutate(
    name = `Compound name`,
    biologicalsource = paste(ifelse(is.na(Genus),
      "",
      Genus
    ),
    ifelse(is.na(Species),
      "",
      Species
    ),
    sep = " "
    ),
    inchi = Inchl,
    smiles = `SMILES (canonical or isomeric)`
  )

data_manipulated$biologicalsource <-
  y_as_na(data_manipulated$biologicalsource, " ")

data_selected <- data_manipulated %>%
  mutate(reference = `Reference_Text Title; Journal; Vol,; Issue; pages; year; type; DOI; author1; authors2; etc.`) %>%
  cSplit("reference", sep = ";") %>%
  select(
    name,
    biologicalsource,
    inchi,
    smiles,
    reference_doi = `DOI   `,
    reference_title = reference_01,
    reference_journal = reference_02
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "cya_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c("reference_doi", "reference_title", "reference_journal")
  )

# exporting
database$writeInterim(data_standard)
