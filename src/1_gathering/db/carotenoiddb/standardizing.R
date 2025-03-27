# title: "CAROTENOIDDB cleaneR"

# loading paths
source("paths.R")
source("r/standardizing_original.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)

# get paths
database <- databases$get("carotenoiddb")

## files
data_original <-
  readr::read_delim(file = gzfile(description = database$sourceFiles$tsv)) |>
  dplyr::mutate_all(as.character)

# manipulating
data_manipulated <- data_original |>
  dplyr::select(
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
    pattern = "(\\))([A-Z])",
    replacement = ")SPLIT\\2",
    x = data_manipulated$biologicalsource
  )

data_manipulated <- data_manipulated |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = "SPLIT",
    stripWhite = FALSE,
    fixed = FALSE
  ) |>
  dplyr::mutate(dplyr::across(.cols = dplyr::everything(), as.character))

data_manipulated_long <- data_manipulated |>
  tidyr::gather(
    8:ncol(data_manipulated),
    key = "n",
    value = "biologicalsource"
  ) |>
  dplyr::group_by(uniqueid) |>
  dplyr::select(-n) |>
  dplyr::distinct(biologicalsource, .keep_all = TRUE) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(biologicalsource) | !n > 1) |>
  dplyr::select(-n) |>
  dplyr::arrange(uniqueid)

data_manipulated_long$reference <-
  gsub(
    pattern = "(Ref.\\d*: )",
    replacement = "SPLIT\\1",
    x = data_manipulated_long$reference
  )

data_manipulated_long_ref <- data_manipulated_long %>%
  splitstackshape::cSplit(
    "reference",
    sep = "SPLIT",
    stripWhite = FALSE,
    fixed = FALSE
  ) %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::gather(8:ncol(.), key = "n", value = "reference") %>%
  dplyr::select(-n) %>%
  group_by(uniqueid) %>%
  dplyr::distinct(biologicalsource, reference, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueid, biologicalsource) %>%
  dplyr::add_count() %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(reference) | !n > 1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(uniqueid)

data_manipulated_long_ref_unique <- data_manipulated_long_ref |>
  dplyr::mutate(
    refnum = stringr::str_extract(
      string = data_manipulated_long_ref$reference,
      pattern = "(Ref.\\d*)"
    ),
    biorefnum = stringr::str_extract(
      string = data_manipulated_long_ref$biologicalsource,
      pattern = "(Ref.\\d*)"
    )
  ) |>
  dplyr::filter(
    refnum == biorefnum |
      is.na(refnum)
  ) |>
  dplyr::mutate(
    reference_unsplittable = gsub(
      pattern = "(Ref.\\d* : )",
      replacement = "",
      x = reference
    )
  ) |>
  dplyr::mutate(
    reference_title = gsub(
      pattern = "\"",
      replacement = "",
      x = stringr::str_extract(
        string = reference_unsplittable,
        pattern = "\".*\""
      )
    ),
    reference_doi_1 = gsub(
      pattern = "\"",
      replacement = "",
      x = stringr::str_extract(
        string = reference_unsplittable,
        pattern = "doi:.*"
      )
    ),
    reference_doi_2 = gsub(
      pattern = "\"",
      replacement = "",
      x = stringr::str_extract(
        string = reference_unsplittable,
        pattern = "DOI:.*"
      )
    )
  ) |>
  dplyr::mutate(
    reference_doi = ifelse(
      test = !is.na(reference_doi_1),
      yes = reference_doi_1,
      no = reference_doi_2
    )
  ) |>
  splitstackshape::cSplit("reference_doi", sep = ",") |>
  splitstackshape::cSplit(
    "biologicalsource",
    sep = " (Ref.",
    fixed = TRUE,
    stripWhite = FALSE
  ) |>
  dplyr::select(
    organism_clean = biologicalsource_1,
    dplyr::everything()
  ) |>
  dplyr::mutate(
    reference_doi = gsub(
      pattern = "doi:",
      replacement = "",
      x = reference_doi_01,
      ignore.case = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_doi = gsub(
      pattern = ". Epub .*",
      replacement = "",
      x = reference_doi,
      ignore.case = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_doi = gsub(
      pattern = " https://doi.org/",
      replacement = "",
      x = reference_doi,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(
    reference_doi = gsub(
      pattern = "http:/​/​dx.​doi.​org/",
      replacement = "",
      x = reference_doi,
      fixed = TRUE
    )
  ) |>
  dplyr::mutate(reference_doi = trimws(reference_doi)) |>
  data.frame() |>
  dplyr::select(
    structure_inchi = inchi,
    structure_smiles = smiles,
    structure_name = name,
    dplyr::everything()
  )

# standardizing
data_standard <-
  standardizing_original(
    data_selected = data_manipulated_long_ref_unique,
    db = "carotenoiddb",
    structure_field = "structure_name",
    organism_field = "organism_clean",
    reference_field = c("reference_title", "reference_doi")
  )

# exporting
database$writeInterim(data_standard)
