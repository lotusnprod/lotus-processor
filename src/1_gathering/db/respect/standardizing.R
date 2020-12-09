# title: "Respect cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(gdata)
library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("respect")

data_standard <- do.call(
  "cbindX",
  lapply(
    database$sourceFiles$tsv,
    function(x) {
      dat <- read.csv(x,
        header = TRUE,
        sep = "\n"
      )
      dat
    }
  )
)

data_transposed <- t(data_standard) %>%
  cSplit(1:ncol(.), sep = ":")

# cleaning
## function
RESPECT_clean <- function(dfsel) {
  df_2 <- dfsel %>%
    select_if(~ sum(!is.na(.)) > 0) %>%
    mutate_all(as.character) %>%
    filter_at(vars(-V1_1), any_vars(. == "SP$SAMPLE")) %>%
    tibble()

  for (i in 1:nrow(df_2)) {
    df_2[i, "biologicalsource_col"] <-
      which(sapply(df_2[i, ], function(x) {
        any(x == "SP$SAMPLE")
      }))
  }

  for (i in 1:nrow(df_2)) {
    df_2[i, "biologicalsource"] <-
      df_2[i, as.numeric((df_2[i, "biologicalsource_col"] + 1))]
  }

  for (i in 1:nrow(df_2)) {
    df_2[i, "inchi_col"] <-
      which(sapply(df_2[i, ], function(x) {
        any(x == "CH$INCHI")
      }))
  }

  for (i in 1:nrow(df_2)) {
    df_2[i, "inchi"] <-
      df_2[i, as.numeric((df_2[i, "inchi_col"] + 1))]
  }

  for (i in 1:nrow(df_2)) {
    df_2[i, "smiles_col"] <-
      which(sapply(df_2[i, ], function(x) {
        any(x == "CH$SMILES")
      }))
  }

  for (i in 1:nrow(df_2)) {
    df_2[i, "smiles"] <-
      df_2[i, as.numeric((df_2[i, "smiles_col"] + 1))]
  }

  df_3 <- df_2

  df_3$biologicalsource <-
    y_as_na(df_3$biologicalsource, "authentic sample")
  df_3$biologicalsource <-
    y_as_na(df_3$biologicalsource, "food sample")
  df_3$biologicalsource <-
    y_as_na(df_3$biologicalsource, "food stuff")
  df_3$biologicalsource <-
    y_as_na(df_3$biologicalsource, "Standard mixture")
  df_3$inchi <- y_as_na(df_3$inchi, "N/A")
  df_3$smiles <- y_as_na(df_3$smiles, "N/A")

  df_4 <- df_3 %>%
    mutate(reference = paste(V3_2, V4_2, sep = "§")) %>%
    select(
      name = 2,
      biologicalsource,
      inchi,
      smiles,
      reference
    )

  df_4$name <- gsub(";.*", "\\1", df_4$name)

  df_4
}

## applying
data_selected <- RESPECT_clean(dfsel = data_transposed)

data_manipulated <- data_selected %>%
  cSplit("reference", sep = "§") %>%
  cSplit("reference_1", sep = ";") %>%
  mutate_all(as.character) %>%
  mutate(
    reference_publishingDetails = paste(reference_1_2,
      reference_1_3,
      reference_1_4,
      sep = ";"
    ),
    reference_pubmed = str_extract(
      string = reference_2,
      pattern = "[0-9]{6,9}"
    )
  ) %>%
  mutate(
    reference_publishingDetails = gsub(
      pattern = ";NA",
      replacement = "",
      x = reference_publishingDetails,
      fixed = TRUE
    )
  ) %>%
  select(
    name,
    biologicalsource,
    inchi,
    smiles,
    reference_pubmed,
    reference_publishingDetails,
    reference_authors = reference_1_1,
    reference_journal = reference_1_2
  ) %>%
  data.frame()

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "res_1",
    structure_field = c("name", "inchi", "smiles"),
    reference_field = c(
      "reference_pubmed",
      "reference_publishingDetails",
      "reference_authors",
      "reference_journal"
    )
  )

# exporting
database$writeInterim(data_standard)