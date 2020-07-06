# title: "SYMMAP cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(readxl)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("symmap")

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

data_original <- fileInZip(inZip = database$sourceFiles$data, varList = "c")

data_bio <- read_excel(database$sourceFiles$bio) %>%
  mutate_all(as.character)

data_chemo <- read_excel(database$sourceFiles$chemo) %>%
  mutate_all(as.character)

# cleaning
data_original$fileName <- gsub('data/', '', data_original$fileName)
data_original$fileName <- gsub('data-', '', data_original$fileName)
data_original$fileName <- gsub('.csv', '', data_original$fileName)
data_original$`Ingredient id` <-
  gsub('SMIT', '', data_original$`Ingredient id`)

data_original <- data_original %>%
  select(
    MOL_id = `Ingredient id`,
    Molecule_name = `Molecule name`,
    Molecule_formula = `Molecule formula`,
    Molecule_weight = `Molecule weight`,
    OB_score = `OB score`,
    Pubchem_id = `PubChem id`,
    CAS_id = `CAS id`,
    Herb_id = fileName
  )

data_chemo$MOL_id <-
  str_pad (data_chemo$MOL_id, width = 5, pad = "0")

colnames(data_bio)[12] <- "TCMID_id_bio"
colnames(data_bio)[13] <- "TCM-ID_id_bio"
colnames(data_bio)[14] <- "TCMSP_id_bio"

data_full <- full_join(data_bio, data_original)

data_full <- full_join(data_chemo, data_full)

# selecting
data_selected <- data_full %>%
  select(
    uniqueid = MOL_id,
    name = Molecule_name,
    pubchem = PubChem_id,
    cas = CAS_id,
    TCMID_id,
    `TCM-ID_id`,
    TCMSP_id,
    TCMID_id_bio,
    `TCM-ID_id_bio`,
    TCMSP_id_bio,
    Latin_name,
    English_name
  ) %>%
  mutate(biologicalsource = paste(Latin_name, English_name, sep = " "))

data_selected$biologicalsource <-
  y_as_na(data_selected$biologicalsource, "NA NA")

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_selected,
    db = "sym_1",
    structure_field = c("name")
  )

# exporting
database$writeInterim(data_standard)
