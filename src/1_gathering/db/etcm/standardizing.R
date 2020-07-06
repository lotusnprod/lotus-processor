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
fileInZip <-
  function(inZip, varList) {
    outFile <- data.frame()
    fileList <- unzip(inZip, list = TRUE)
    for (i in 1:nrow(fileList)) {
      if (grepl("csv", fileList[i, 1])) {
        oFa <- read_delim(
          unz(inZip, fileList[i, 1]),
          delim = ",",
          escape_double = FALSE,
          trim_ws = TRUE
        )
        oFa$fileName <- fileList[i, 1]
        outFile <-
          rbind(outFile, oFa)
      }
    }
    return(outFile)
  }

data_original <-
  fileInZip(inZip = database$sourceFiles$data, varList = "c")

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
  standardizing_original(
    data_selected = data_selected_long,
    db = "etc_1",
    structure_field = "name",
    reference_field = c()
  )

# exporting
database$writeInterim(data_standard)
