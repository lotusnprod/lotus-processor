# title: "ETCM cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("etcm")

## files
data_original <- do.call("rbind",
                         lapply(database$sourceFiles$tsv,
                                function(x) {
                                  dat <- read_delim(
                                    file = x,
                                    delim = ",",
                                    escape_double = FALSE,
                                    trim_ws = TRUE
                                  ) %>%
                                    mutate_all(as.character)
                                  dat$fileName <-
                                    tools::file_path_sans_ext(basename(x))
                                  dat
                                })) %>%
  mutate_all(as.character)

# cleaning
data_original$fileName <-
  gsub('tableExport-', '', data_original$fileName)

data_wide <-
  pivot_wider(data_original, names_from = X1, values_from = X2) %>%
  select(
    herbid = fileName,
    latin = `Herb Name in Ladin`,
    family = `Description in English`,
    name = Components,
    flavor = Flavor
  )

data_wide$family <- y_as_na(x = data_wide$family , y = "")

data_wide <- data_wide %>%
  filter(!is.na(family))

data_wide_2 <- data_wide %>%
  cSplit("name", sep = ", ") %>%
  pivot_longer(5:ncol(.))

# selecting
data_selected <- data_wide_2 %>%
  select(herbid,
         latin,
         family,
         flavor,
         name = value) %>%
  filter(!is.na(name))

data_selected_long <- data_selected %>%
  separate(
    col = "name",
    into = c("name01", "name02", "name03", "name04", "name05", "name06"),
    sep = ",(?=[A-Z])"
  ) %>%
  pivot_longer(5:ncol(.)) %>%
  filter(!is.na(value)) %>%
  select(
    herbid,
    name_latin = latin,
    name_family = family,
    flavor,
    name = value
  ) %>%
  mutate(biologicalsource = paste(name_latin[!is.na(name_latin)],
                                  name_family[!is.na(name_family)],
                                  sep = " ")) %>%
  select(name,
         biologicalsource)

# standardizing
data_standard <-
  standardizing_original(data_selected = data_selected_long,
                         db = "etc_1",
                         structure_field = "name")

# exporting
database$writeInterim(data_standard)
