# title: "ALKAMID cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(splitstackshape)
library(readr)
library(tidyr)

# get paths
database <- databases$get("alkamid")

## db
data_original <- read_delim(
  file = gzfile(database$sourceFiles$tsv)
) %>%
  mutate_all(as.character) %>%
  mutate(id = row.names(.))

## ref
ref_original <- read_delim(
  file = gzfile(database$sourceFiles$tsvRef)
) %>%
  mutate_all(as.character)

# cleaning
cleaning_alkamid <- function(x) {
  df1 <- x %>%
    select(32:ncol(.))

  df11 <- df1 %>%
    filter(V1_032 == "No .hin file found to display.")

  df12 <- df1 %>%
    filter(V1_032 != "No .hin file found to display.")

  for (i in (ncol(df11):9)) {
    df11[, i] <- df11[, i - 8]
  }

  df2 <- full_join(df12, df11) %>%
    select(11:ncol(.))

  df21 <- df2 %>%
    filter(!grepl("Alkamid", V1_042))

  df22 <- df2 %>%
    filter(grepl("Alkamid", V1_042))

  for (i in (ncol(df22):2)) {
    df22[, i] <- df22[, i - 1]
  }

  df3 <- full_join(df22, df21) %>%
    select(
      19,
      25:ncol(.)
    )

  df3$biologicalsource <-
    apply(df3[, 5:ncol(df3)], 1, paste, collapse = "|")

  df3$biologicalsource <-
    gsub("Functionality.*", "", df3$biologicalsource, fixed = FALSE)

  df4 <- df3 %>%
    select(
      name = V1_060,
      biologicalsource
    ) %>%
    filter(!grepl("NA", biologicalsource)) %>%
    filter(!grepl("MW", biologicalsource)) %>%
    cSplit("biologicalsource", sep = "|")

  df41 <- df4 %>%
    filter(!grepl("Tribe", biologicalsource_01)) %>%
    tibble()

  df42 <- df4 %>%
    filter(grepl("Tribe", biologicalsource_01))

  for (i in (2:(ncol(df41) - 1))) {
    df41[, i] <- df41[, i + 1]
  }

  df5 <- full_join(df42, df41) %>%
    pivot_longer(
      2:ncol(.),
      names_to = c(".value", "level"),
      names_sep = "_",
      values_to = "newcanonicalname",
      values_drop_na = TRUE
    ) %>%
    mutate(value = as.character(biologicalsource))

  shift <- function(x, n) {
    c(x[-(seq(n))], rep(NA, n))
  }

  df5$value <- shift(df5$value, 1)

  df6 <- df5 %>%
    filter(row_number() %% 2 == 1) %>%
    select(-level) %>%
    pivot_wider(names_from = biologicalsource, values_from = value) %>%
    unnest(Genus) %>%
    unnest(Species) %>%
    mutate(biologicalsource = paste(Genus, Species, sep = " ")) %>%
    select(
      name,
      biologicalsource
    )

  df7 <- df6 %>%
    cSplit("name", sep = "Emperical formula", stripWhite = FALSE) %>%
    cSplit("name_1", sep = "SMILES", stripWhite = FALSE) %>%
    cSplit("name_1_1", sep = "Trivial name", stripWhite = FALSE) %>%
    cSplit("name_1_1_1", sep = "IUPAC name", stripWhite = FALSE) %>%
    cSplit("name_1_1_1_1", sep = "Chemical name", stripWhite = FALSE) %>%
    cSplit(
      "name_1_1_1_1_1",
      sep = "	Alkamid moleculeID",
      stripWhite = FALSE
    ) %>%
    select(
      uniqueid = name_1_1_1_1_1_1,
      name = name_1_1_2,
      smiles = name_1_2,
      organism_clean = biologicalsource
    )
  return(df7)
}

# selecting
data_selected <- cleaning_alkamid(x = data_original)

# making names compatible
ref_prepared <- ref_original %>%
  mutate(uniqueid = paste0("Alkamid moleculeID", entry_id)) %>%
  select(
    uniqueid,
    reference_authors = `Author(s)`,
    reference_journal = Journal,
    reference_title = Title
  )

# joining references
data_referenced <- left_join(data_selected, ref_prepared) %>%
  group_by(uniqueid, organism_clean) %>%
  select(
    structure_smiles = smiles,
    structure_name = name,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_referenced,
    db = "alkamid",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_authors",
      "reference_journal",
      "reference_title"
    )
  )

# exporting
database$writeInterim(data_standard)
