# title: "TMDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

groundhog.library(splitstackshape, date = groundhog.day)
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(vroom, date = groundhog.day)

# get paths
database <- databases$get("tmdb")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  quote = ""
)

# pivoting
data_pivoted <- data_original %>%
  mutate(level = as.numeric(gl(nrow(.) / 28, 28))) %>%
  group_by(level) %>%
  pivot_wider(names_from = 1, values_from = 2) %>%
  unnest() %>%
  ungroup()

# selecting
data_selected <- data_pivoted %>%
  select(
    name = `Entry name`,
    biologicalsource = `Latin name`,
    reference_publishingDetails = References
  ) %>%
  cSplit("reference_publishingDetails",
    sep = ";",
    direction = "long"
  ) %>%
  mutate(
    reference_publishingDetails = gsub(
      pattern = "\\(.*\\D.*\\)",
      replacement = "",
      x = reference_publishingDetails
    )
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "tmd_1",
    structure_field = c("name"),
    reference_field = c("reference_publishingDetails")
  )

data_standard[] <-
  lapply(data_standard, function(x) {
    gsub("Not Available", NA, x)
  })

# exporting
database$writeInterim(data_standard)
