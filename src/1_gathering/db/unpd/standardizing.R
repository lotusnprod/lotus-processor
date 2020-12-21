# title: "UNPD cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)

# get paths
database <- databases$get("unpd")

## files
data_original <- read_delim(
  file = gzfile(pathDataExternalDbSourceUnpdIntegrated),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

# selecting
## atomizing references
data_selected <- data_original %>%
  mutate(reference = gsub("(\\(\\d+).\\s", "|", ref)) %>%
  cSplit("reference", sep = "|", direction = "long") %>%
  mutate_all(as.character) %>%
  cSplit(
    "reference",
    sep = "; ; ; ",
    stripWhite = FALSE,
    fixed = TRUE,
    direction = "long"
  ) %>%
  mutate_all(as.character) %>%
  select(
    biologicalsource = ln_reduced,
    reference_original = reference,
    inchi = InChI,
    smiles = SMILES
  ) %>%
  data.frame()

data_manipulated <- data_selected %>%
  mutate(
    reference_external = ifelse(
      test = reference_original == "Retrieved from CNPD",
      yes = reference_original,
      no = NA
    ),
    reference_original = gsub(
      pattern = "Retrieved from CNPD",
      replacement = "",
      x = reference_original,
      fixed = TRUE
    )
  ) %>%
  mutate(reference_original = y_as_na(reference_original, y = "")) %>%
  mutate(reference_original = ifelse(
    test = !grepl(pattern = "[^ -~]", x = reference_original),
    yes = reference_original,
    no = NA
  )) %>%
  mutate(
    reference_split =
      ifelse(
        test = grepl(
          pattern = ".*et al.",
          x = reference_original
        ),
        yes = trimws(x = sub(
          pattern = "^.;",
          replacement = "",
          x = sub(
            pattern = "^ \\\\",
            replacement = "",
            x = sub(
              pattern = "^\\\\",
              replacement = "",
              x = sub(
                pattern = "^ \\.",
                replacement = "",
                x = sub(
                  pattern = "^;",
                  replacement = "",
                  x = sub(
                    pattern = ".*et al.",
                    replacement = "",
                    x = reference_original
                  )
                )
              )
            )
          )
        )),
        no = NA
      )
  ) %>%
  data.frame()

# reverse_words <- function(string) {
#   # split string by blank spaces
#   string_split = strsplit(as.character(string), split = " ")
#   # how many split terms?
#   string_length = length(string_split[[1]])
#   # decide what to do
#   if (string_length == 1) {
#     # one word (do nothing)
#     reversed_string = string_split[[1]]
#   } else {
#     # more than one word (collapse them)
#     reversed_split = string_split[[1]][string_length:1]
#     reversed_string = paste(reversed_split, collapse = " ")
#   }
#   # output
#   return(reversed_string)
# }

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated,
    db = "unp_1",
    structure_field = c("inchi", "name", "smiles"),
    reference_field = c(
      "reference_original",
      "reference_external",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)
