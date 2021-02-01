# title: "CAROTENOIDDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(splitstackshape)
library(tidyverse)
library(vroom)

# get paths
database <- databases$get("carotenoiddb")

## files
data_original <- vroom(
  file = gzfile(database$sourceFiles$tsv),
  delim = "\t",
  quote = ""
) %>%
  mutate_all(as.character)

# manipulating
data_manipulated <- data_original %>%
  select(
    uniqueid = Entry,
    name = `Carotenoid Name`,
    inchi = InChI,
    inchikey = InChIKey,
    smiles = `Canonical SMILES`,
    cas = CAS,
    biologicalsource = `Source organisms`,
    reference = References
  )

data_manipulated$biologicalsource <-
  gsub(
    "(\\))([A-Z])",
    ")SPLIT\\2",
    data_manipulated$biologicalsource
  )

data_manipulated <- data_manipulated %>%
  cSplit(
    "biologicalsource",
    sep = "SPLIT",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  mutate(across(.cols = everything(), as.character))

data_manipulated_long <- data_manipulated %>%
  gather(8:ncol(.),
    key = "n",
    value = "biologicalsource"
  ) %>%
  group_by(uniqueid) %>%
  select(-n) %>%
  distinct(biologicalsource, .keep_all = TRUE) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(biologicalsource) | !n > 1) %>%
  select(-n) %>%
  arrange(uniqueid)

data_manipulated_long$reference <-
  gsub("(Ref.\\d* : )", "SPLIT\\1", data_manipulated_long$reference)

data_manipulated_long_ref <- data_manipulated_long %>%
  cSplit("reference",
    sep = "SPLIT",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  mutate_all(as.character) %>%
  gather(8:ncol(.),
    key = "n",
    value = "reference"
  ) %>%
  select(-n) %>%
  group_by(uniqueid) %>%
  distinct(biologicalsource, reference, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(uniqueid, biologicalsource) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(reference) | !n > 1) %>%
  select(-n) %>%
  arrange(uniqueid)

data_manipulated_long_ref_unique <- data_manipulated_long_ref %>%
  mutate(
    refnum = str_extract(data_manipulated_long_ref$reference, "(Ref.\\d*)"),
    biorefnum = str_extract(data_manipulated_long_ref$biologicalsource, "(Ref.\\d*)")
  ) %>%
  filter(refnum == biorefnum |
    is.na(refnum)) %>%
  mutate(reference_unsplittable = gsub(
    pattern = "(Ref.\\d* : )",
    replacement = "",
    x = reference
  )) %>%
  mutate(
    reference_title = gsub(
      pattern = "\"",
      replacement = "",
      x = str_extract(
        string = reference_unsplittable,
        pattern = "\".*\""
      )
    ),
    reference_doi_1 =
      gsub(
        pattern = "\"",
        replacement = "",
        x = str_extract(
          string = reference_unsplittable,
          pattern = "doi:.*"
        )
      ),
    reference_doi_2 =
      gsub(
        pattern = "\"",
        replacement = "",
        x = str_extract(
          string = reference_unsplittable,
          pattern = "DOI:.*"
        )
      )
  ) %>%
  mutate(reference_doi = ifelse(
    test = !is.na(reference_doi_1),
    yes = reference_doi_1,
    no = reference_doi_2
  )) %>%
  cSplit("reference_doi", sep = ",") %>%
  cSplit(
    "biologicalsource",
    sep = " (Ref.",
    fixed = TRUE,
    stripWhite = FALSE
  ) %>%
  select(
    organism_clean = biologicalsource_1,
    everything()
  ) %>%
  mutate(
    reference_doi = gsub(
      pattern = "doi:",
      replacement = "",
      x = reference_doi_01,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    reference_doi = gsub(
      pattern = ". Epub .*",
      replacement = "",
      x = reference_doi,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    reference_doi = gsub(
      pattern = " https://doi.org/",
      replacement = "",
      x = reference_doi,
      fixed = TRUE
    )
  ) %>%
  mutate(
    reference_doi = gsub(
      pattern = "http:/​/​dx.​doi.​org/",
      replacement = "",
      x = reference_doi,
      fixed = TRUE
    )
  ) %>%
  mutate(reference_doi = trimws(reference_doi)) %>%
  data.frame() %>%
  select(
    structure_inchi = inchi,
    structure_smiles = smiles,
    structure_name = name,
    everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated_long_ref_unique,
    db = "car_1",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c("reference_title", "reference_doi")
  )

# exporting
database$writeInterim(data_standard)
