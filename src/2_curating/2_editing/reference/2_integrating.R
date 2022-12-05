source("r/log_debug.R")
log_debug("This script integrates all reference translations together")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

source("r/split_data_table_quote.R")

log_debug("... files ...")
log_debug("... DOI")
dataDoi <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferenceDoi,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::select(
    referenceOriginal = referenceOriginal_doi,
    doi_doi = referenceTranslatedDoi,
    journal_doi = referenceTranslatedJournal,
    title_doi = referenceTranslatedTitle,
    date_doi = referenceTranslatedDate,
    author_doi = referenceTranslatedAuthor,
    scoreCrossref_doi = referenceTranslationScoreCrossref,
    scoreDistance_doi = referenceTranslationScoreDistance
  ) |>
  dplyr::mutate_all(as.character) |>
  tidyr::pivot_longer(
    cols = 2:8,
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) |>
  dplyr::mutate(level = 1) |>
  dplyr::mutate_all(as.character)

log_debug("... original references")
dataOriginal <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferenceOriginal,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  dplyr::select(
    referenceOriginal = referenceOriginal_original,
    doi_original = referenceTranslatedDoi,
    journal_original = referenceTranslatedJournal,
    title_original = referenceTranslatedTitle,
    date_original = referenceTranslatedDate,
    author_original = referenceTranslatedAuthor,
    scoreCrossref_original = referenceTranslationScoreCrossref,
    scoreDistance_original = referenceTranslationScoreDistance
  ) %>%
  dplyr::group_by(referenceOriginal) %>%
  dplyr::mutate(level = dplyr::row_number()) %>%
  dplyr::relocate(level, .after = referenceOriginal) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  dplyr::mutate_all(as.character)

log_debug("... PMID")
dataPubmed <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferencePubmed,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  select(
    referenceOriginal = referenceOriginal_pubmed,
    doi_pubmed = referenceTranslatedDoi,
    journal_pubmed = referenceTranslatedJournal,
    title_pubmed = referenceTranslatedTitle,
    date_pubmed = referenceTranslatedDate,
    author_pubmed = referenceTranslatedAuthor,
    scoreCrossref_pubmed = referenceTranslationScoreCrossref,
    scoreDistance_pubmed = referenceTranslationScoreDistance
  ) |>
  dplyr::mutate_all(as.character) |>
  tidyr::pivot_longer(
    cols = 2:8,
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) |>
  dplyr::mutate(level = 1) |>
  dplyr::mutate_all(as.character)

log_debug("... publishing details")
dataPublishingDetails <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferencePublishingDetails,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  dplyr::select(
    referenceOriginal = referenceOriginal_publishingDetails,
    doi_publishingDetails = referenceTranslatedDoi,
    journal_publishingDetails = referenceTranslatedJournal,
    title_publishingDetails = referenceTranslatedTitle,
    date_publishingDetails = referenceTranslatedDate,
    author_publishingDetails = referenceTranslatedAuthor,
    scoreCrossref_publishingDetails = referenceTranslationScoreCrossref,
    scoreDistance_publishingDetails = referenceTranslationScoreDistance
  ) %>%
  dplyr::group_by(referenceOriginal) %>%
  dplyr::mutate(level = dplyr::row_number()) %>%
  dplyr::relocate(level, .after = referenceOriginal) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  dplyr::mutate_all(as.character)

log_debug("... titles")
dataTitle <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferenceTitle,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  dplyr::select(
    referenceOriginal = referenceOriginal_title,
    doi_title = referenceTranslatedDoi,
    journal_title = referenceTranslatedJournal,
    title_title = referenceTranslatedTitle,
    date_title = referenceTranslatedDate,
    author_title = referenceTranslatedAuthor,
    scoreCrossref_title = referenceTranslationScoreCrossref,
    scoreDistance_title = referenceTranslationScoreDistance
  ) %>%
  dplyr::group_by(referenceOriginal) %>%
  dplyr::mutate(level = dplyr::row_number()) %>%
  dplyr::relocate(level, .after = referenceOriginal) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  dplyr::mutate_all(as.character)

log_debug("... split")
dataSplit <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedReferenceSplit,
    delim = "\t",
    col_types = cols(.default = "c")
  ) %>%
  dplyr::select(
    referenceOriginal = referenceOriginal_split,
    doi_split = referenceTranslatedDoi,
    journal_split = referenceTranslatedJournal,
    title_split = referenceTranslatedTitle,
    date_split = referenceTranslatedDate,
    author_split = referenceTranslatedAuthor,
    scoreCrossref_split = referenceTranslationScoreCrossref,
    scoreDistance_split = referenceTranslationScoreDistance
  ) %>%
  dplyr::group_by(referenceOriginal) %>%
  dplyr::mutate(level = dplyr::row_number()) %>%
  dplyr::relocate(level, .after = referenceOriginal) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_all(as.character) %>%
  tidyr::pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  dplyr::mutate_all(as.character)

log_debug("... full references")
dataFull <-
  readr::read_delim(
    file = pathDataInterimTablesOriginalReferenceFull,
    delim = "\t",
    col_types = cols(.default = "c")
  )

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("...  cleaned organisms")
  dataCleanedOrganismManipulated_old <-
    readr::read_delim(
      file = pathDataInterimDictionariesOrganismDictionary,
      delim = "\t",
      col_types = cols(.default = "c")
    ) |>
    dplyr::mutate(
      organismDetected =
        word(organismDetected, 1, 2)
    ) |>
    dplyr::distinct(
      organismType,
      organismValue,
      organismDetected
    )

  dataCleanedOrganismManipulated_new <-
    readr::read_delim(
      file = pathDataInterimTablesProcessedOrganismFinal,
      delim = "\t",
      col_types = cols(.default = "c")
    ) |>
    dplyr::mutate(
      organismDetected =
        word(organismDetected, 1, 2)
    ) |>
    dplyr::distinct(
      organismType,
      organismValue,
      organismDetected
    )

  dataCleanedOrganismManipulated <-
    dplyr::bind_rows(
      dataCleanedOrganismManipulated_old,
      dataCleanedOrganismManipulated_new
    ) |>
    dplyr::distinct()
}

if (!file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("... cleaned organisms")
  dataCleanedOrganismManipulated <-
    readr::read_delim(
      file = pathDataInterimTablesProcessedOrganismFinal,
      delim = "\t",
      col_types = cols(.default = "c")
    ) |>
    dplyr::mutate(
      organismDetected =
        stringr::word(
          string = organismDetected,
          start = 1,
          end = 2
        )
    ) |>
    dplyr::distinct(
      organismType,
      organismValue,
      organismDetected
    )
}

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  log_debug("... reference dictionary, this may take a while")
  referenceDictionary <-
    readr::read_delim(
      file = pathDataInterimDictionariesReferenceDictionary,
      delim = "\t",
      col_types = cols(.default = "c")
    )
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  log_debug("... reference organism dictionary")
  referenceOrganismDictionary <-
    readr::read_delim(
      file = pathDataInterimDictionariesReferenceOrganismDictionary,
      delim = "\t",
      col_types = cols(.default = "c")
    )
}

log_debug("joining ...")
log_debug("... all reference types")
dataCrossref <- dplyr::bind_rows(
  dataDoi,
  dataOriginal,
  dataPublishingDetails,
  dataPubmed,
  dataSplit,
  dataTitle
) |>
  dplyr::filter(!is.na(referenceOriginal)) |>
  dplyr::filter(!is.na(referenceTranslatedValue)) |>
  dplyr::distinct(
    referenceOriginal,
    referenceTranslatedType,
    origin,
    referenceTranslatedValue,
    level
  )

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  dataCrossref <- dplyr::bind_rows(dataCrossref, referenceDictionary)
}

dataCrossref <- dataCrossref |>
  dplyr::filter(!is.na(referenceOriginal)) |>
  dplyr::filter(!is.na(referenceTranslatedValue)) |>
  dplyr::distinct(
    referenceOriginal,
    referenceTranslatedType,
    origin,
    referenceTranslatedValue,
    level
  )

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  dataFull <- dataFull |>
    dplyr::anti_join(referenceOrganismDictionary)
}

log_debug("... with organisms")
dataTranslated <-
  dplyr::left_join(dataFull, dataCleanedOrganismManipulated) |>
  dplyr::filter(!is.na(referenceValue)) |>
  dplyr::distinct(
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    organismDetected
  ) |>
  data.table::data.table()

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimDictionariesReference)

log_debug("exporting, this may take a while if running full mode")
log_debug(pathDataInterimTablesTranslatedReference)
split_data_table_quote(
  x = dataTranslated,
  no_rows_per_frame = 100000,
  text = "",
  path_to_store = pathDataInterimTablesTranslatedReference
)

log_debug(pathDataInterimDictionariesReferenceDictionary)
readr::write_delim(
  x = dataCrossref,
  delim = "\t",
  file = pathDataInterimDictionariesReferenceDictionary,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
