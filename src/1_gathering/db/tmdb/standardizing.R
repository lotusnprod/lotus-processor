# title: "TMDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("tmdb")

## files
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# pivoting
data_pivoted <- data_original %>%
  mutate(level = as.numeric(gl(nrow(.) / 28, 28))) %>%
  group_by(level) %>%
  pivot_wider(names_from = 1, values_from = 2) %>%
  unnest() %>%
  ungroup()


# selecting
data_selected <- data_pivoted %>%
  select(name = `Entry name`,
         biologicalsource = `Latin name`,
         reference = References)


# standardizing
data_standard <-
  standardizing_original(data_selected = data_selected,
                         db = "tmd_1",
                         structure_field = "name")

data_standard[] <-
  lapply(data_standard, function(x)
    gsub("Not Available", NA, x))

# exporting
database$writeInterim(data_standard)