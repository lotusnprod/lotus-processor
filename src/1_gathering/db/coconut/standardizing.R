# title: "COCONUT cleaneR"

# loading paths
source("paths.R")
source("r/y_as_na.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("coconut")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv))

"%ni%" <- Negate("%in%")

# selecting
data_selected <- data_original |>
  dplyr::select(
    name,
    inchi = inchi,
    smiles = SMILES,
    biologicalsource = textTaxa,
    reference = citationDOI,
    reference_external = found_in_databases
  ) %>%
  dplyr::filter(biologicalsource != "[notax]") %>%
  dplyr::filter(biologicalsource != "[]") %>%
  dplyr::filter(reference != "[]") %>%
  dplyr::mutate(
    biologicalsource = gsub(
      pattern = "\\[",
      replacement = "",
      x = biologicalsource
    ),
    biologicalsource = gsub(
      pattern = "\\]",
      replacement = "",
      x = biologicalsource
    ),
    biologicalsource = gsub(
      pattern = "\"",
      replacement = "",
      x = biologicalsource
    ),
    reference = gsub(
      pattern = "\\[",
      replacement = "",
      x = reference
    ),
    reference = gsub(
      pattern = "\\]",
      replacement = "",
      x = reference
    ),
    reference = gsub(
      pattern = "\"",
      replacement = "",
      x = reference
    ),
    # reference_external = gsub(pattern = "\\[", replacement = "", x = reference_external),
    # reference_external = gsub(pattern = "\\]", replacement = "", x = reference_external),
    # reference_external = gsub(pattern = "\"", replacement = "", x = reference_external)
  )

data_counted <- data_selected |>
  dplyr::rowwise() |>
  dplyr::mutate(
    nrefs = stringr::str_count(string = reference, pattern = ", "),
    norganisms = stringr::str_count(string = biologicalsource, pattern = ", ")
  )

## right organism, right ref
data_nosplit <- data_counted |>
  dplyr::filter(nrefs == 0)

## right organism, right ref
data_easy <- data_counted |>
  dplyr::filter(nrefs != 0 & nrefs == norganisms)

## combination of all organisms vs all references since we do not know which one is good
data_hard <- data_counted |>
  dplyr::filter(nrefs != 0 & nrefs != norganisms)

data_hard_harborne <- data_hard |>
  dplyr::filter(
    grepl(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      x = reference
    )
  )

data_hard_harborne_treated <- data_hard_harborne |>
  dplyr::mutate(
    reference = gsub(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      replacement = "",
      x = reference
    )
  ) |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::group_by(smiles) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) |>
  dplyr::select(-n) |>
  splitstackshape::cSplit("reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::filter(!grepl(pattern = "^,", x = reference)) |>
  dplyr::filter(!grepl(pattern = "^\\.", x = reference)) |>
  dplyr::filter(!grepl(pattern = "^[0-9]{1-3},", x = reference)) |>
  dplyr::filter(!grepl(pattern = "^\\([0-9]{4}\\)", x = reference)) |>
  dplyr::filter(stringr::str_length(string = reference) >= 8) |>
  dplyr::add_count(reference) |>
  dplyr::filter(
    n <= 5400 |
      grepl(
        pattern = "Molecules_2015;20(8):15330-42",
        x = reference
      ) |
      grepl(
        pattern = "Lu,Phytochem.,59,(2002),117",
        x = reference
      )
  ) |>
  dplyr::select(-nrefs, -norganisms, -n) |>
  data.frame()

data_hard_2 <- data_hard |>
  dplyr::filter(
    !grepl(
      pattern = "Harborne, The Handbook of Natural Flavonoids, [0-9], \\([0-9]{4}\\), [0-9]{1-3}",
      x = reference
    )
  )

data_nosplit_treated <- data_nosplit |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::group_by(smiles) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) |>
  dplyr::select(-n) |>
  splitstackshape::cSplit("reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::filter(stringr::str_length(string = reference) >= 8) |>
  dplyr::select(-nrefs, -norganisms) |>
  data.frame()

data_easy_treated <- data_easy |>
  splitstackshape::cSplit("biologicalsource",
    sep = ", ",
    fixed = TRUE
  ) %>%
  splitstackshape::cSplit("reference",
    sep = ", ",
    fixed = TRUE
  ) %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::pivot_longer(
    cols = 7:(ncol(.)),
    names_to = c("type", "number"),
    names_sep = "_"
  ) |>
  dplyr::filter(!is.na(value)) %>%
  tidyr::pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::group_by(inchi) %>%
  tidyr::fill(reference, .direction = "downup") %>%
  dplyr::group_by(inchi) %>%
  dplyr::add_count() %>%
  dplyr::ungroup() %>%
  dplyr::filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) %>%
  dplyr::select(-n) %>%
  dplyr::filter(stringr::str_length(string = reference) >= 8) %>%
  dplyr::select(-nrefs, -norganisms, -number) %>%
  data.frame()

data_hard_treated <- data_hard_2 %>%
  splitstackshape::cSplit(
    "biologicalsource",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) %>%
  dplyr::group_by(smiles) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(n == 1 |
    biologicalsource %ni% c("plants", "marine", "fungi")) |>
  dplyr::select(-n) |>
  splitstackshape::cSplit("reference",
    sep = ", ",
    direction = "long",
    fixed = TRUE
  ) |>
  dplyr::filter(stringr::str_length(string = reference) >= 8) |>
  dplyr::add_count(reference) |>
  dplyr::filter(n <= 15000) |>
  dplyr::select(-nrefs, -norganisms, -n) |>
  data.frame()

data_treated <- rbind(
  data_nosplit_treated,
  data_easy_treated,
  data_hard_harborne_treated,
  data_hard_treated
) |>
  dplyr::filter(grepl(pattern = "[0-9]", x = reference))

data_corrected <- data_treated |>
  dplyr::mutate(
    reference_split = reference,
    reference_authors = reference
  ) |>
  dplyr::mutate(n = nchar(reference)) |>
  dplyr::mutate(
    reference_pubmed = ifelse(
      test = n == 8,
      yes = stringr::str_extract(string = reference, pattern = "[0-9]{8}"),
      no = NA
    ),
    reference_doi = ifelse(
      test = n >= 9,
      yes = stringr::str_extract(string = reference, pattern = "^10.*"),
      no = NA
    ),
    reference_split = ifelse(
      test = n >= 20,
      yes = stringr::str_extract(string = reference_split, pattern = "^[A-Z].*"),
      no = NA
    ),
    reference_authors = ifelse(
      test = n < 20,
      yes = stringr::str_extract(string = reference_authors, pattern = "^[A-Z].*"),
      no = NA
    )
  ) |>
  dplyr::select(-n)

data_corrected$name <- y_as_na(data_corrected$name, "")
data_corrected$inchi <- y_as_na(data_corrected$inchi, "")
data_corrected$biologicalsource <-
  y_as_na(data_corrected$biologicalsource, "")
data_corrected$reference_external <-
  y_as_na(data_corrected$reference_external, "")

data_corrected$name <- y_as_na(data_corrected$name, "NA")
data_corrected$inchi <- y_as_na(data_corrected$inchi, "NA")
data_corrected$biologicalsource <-
  y_as_na(data_corrected$biologicalsource, "NA")
data_corrected$reference_authors <-
  y_as_na(data_corrected$reference_authors, "NA")
data_corrected$reference_external <-
  y_as_na(data_corrected$reference_external, "NA")

data_corrected$biologicalsource <-
  gsub(
    pattern = "_",
    replacement = " ",
    x = data_corrected$biologicalsource,
    fixed = TRUE
  )

findCapitals_1 <- data_corrected$biologicalsource %>%
  data.frame() %>%
  dplyr::filter(!is.na(.)) %>%
  dplyr::mutate(nchar = nchar(.)) %>%
  dplyr::mutate(first = substr(., start = 1, stop = 1)) %>%
  dplyr::mutate(last = substr(., start = nchar, stop = nchar)) %>%
  dplyr::select(
    x = ".",
    nchar,
    first,
    last
  ) %>%
  dplyr::distinct()

findCapitals_2 <- data_corrected$biologicalsource %>%
  data.frame() %>%
  dplyr::filter(!is.na(.)) %>%
  dplyr::mutate(nchar = nchar(.)) %>%
  dplyr::mutate(first = substr(., start = 1, stop = 1)) %>%
  dplyr::mutate(last = substr(., start = nchar, stop = nchar)) %>%
  dplyr::select(
    y = ".",
    nchar,
    first,
    last
  ) %>%
  dplyr::distinct()

findCapitals_3 <-
  dplyr::full_join(findCapitals_1, findCapitals_2) |>
  dplyr::filter(tolower(x) == tolower(y)) |>
  dplyr::filter(x != y) |>
  dplyr::mutate(
    capitals_x = stringr::str_count(string = x, pattern = "[A-Z]"),
    capitals_y = stringr::str_count(string = y, pattern = "[A-Z]")
  ) |>
  dplyr::filter(capitals_x > capitals_y)

data_corrected_capitals <- dplyr::left_join(data_corrected,
  findCapitals_3,
  by = c("biologicalsource" = "x")
) |>
  dplyr::mutate(organism_clean = biologicalsource) |>
  dplyr::select(-y, -nchar, -first, -last, -capitals_x, -capitals_y) |>
  dplyr::select(
    structure_inchi = inchi,
    structure_smiles = smiles,
    structure_name = name,
    dplyr::everything()
  ) |>
  dplyr::filter(!grepl(
    pattern = "$",
    x = organism_clean,
    fixed = TRUE
  ))

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_corrected_capitals,
    db = "coconut",
    structure_field = "structure_smiles",
    organism_field = "organism_clean",
    reference_field = c(
      "reference_doi",
      "reference_pubmed",
      "reference_authors",
      "reference_split"
    )
  )

# exporting
database$writeInterim(data_standard)
