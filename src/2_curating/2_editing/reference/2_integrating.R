source("r/log_debug.R")
log_debug("This script integrates all reference translations together")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

log_debug("... files ...")
log_debug("... DOI")
dataDoi <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferenceDoi,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_doi,
    doi_doi = referenceTranslatedDoi,
    journal_doi = referenceTranslatedJournal,
    title_doi = referenceTranslatedTitle,
    date_doi = referenceTranslatedDate,
    author_doi = referenceTranslatedAuthor,
    scoreCrossref_doi = referenceTranslationScoreCrossref,
    scoreDistance_doi = referenceTranslationScoreDistance
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate(level = 1) %>%
  mutate_all(as.character)

log_debug("... original references")
dataOriginal <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferenceOriginal,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_original,
    doi_original = referenceTranslatedDoi,
    journal_original = referenceTranslatedJournal,
    title_original = referenceTranslatedTitle,
    date_original = referenceTranslatedDate,
    author_original = referenceTranslatedAuthor,
    scoreCrossref_original = referenceTranslationScoreCrossref,
    scoreDistance_original = referenceTranslationScoreDistance
  ) %>%
  group_by(referenceOriginal) %>%
  mutate(level = row_number()) %>%
  relocate(level, .after = referenceOriginal) %>%
  ungroup() %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate_all(as.character)

log_debug("... PMID")
dataPubmed <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferencePubmed,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_pubmed,
    doi_pubmed = referenceTranslatedDoi,
    journal_pubmed = referenceTranslatedJournal,
    title_pubmed = referenceTranslatedTitle,
    date_pubmed = referenceTranslatedDate,
    author_pubmed = referenceTranslatedAuthor,
    scoreCrossref_pubmed = referenceTranslationScoreCrossref,
    scoreDistance_pubmed = referenceTranslationScoreDistance
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate(level = 1) %>%
  mutate_all(as.character)

log_debug("... publishing details")
dataPublishingDetails <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferencePublishingDetails,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_publishingDetails,
    doi_publishingDetails = referenceTranslatedDoi,
    journal_publishingDetails = referenceTranslatedJournal,
    title_publishingDetails = referenceTranslatedTitle,
    date_publishingDetails = referenceTranslatedDate,
    author_publishingDetails = referenceTranslatedAuthor,
    scoreCrossref_publishingDetails = referenceTranslationScoreCrossref,
    scoreDistance_publishingDetails = referenceTranslationScoreDistance
  ) %>%
  group_by(referenceOriginal) %>%
  mutate(level = row_number()) %>%
  relocate(level, .after = referenceOriginal) %>%
  ungroup() %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate_all(as.character)

log_debug("... titles")
dataTitle <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferenceTitle,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_title,
    doi_title = referenceTranslatedDoi,
    journal_title = referenceTranslatedJournal,
    title_title = referenceTranslatedTitle,
    date_title = referenceTranslatedDate,
    author_title = referenceTranslatedAuthor,
    scoreCrossref_title = referenceTranslationScoreCrossref,
    scoreDistance_title = referenceTranslationScoreDistance
  ) %>%
  group_by(referenceOriginal) %>%
  mutate(level = row_number()) %>%
  relocate(level, .after = referenceOriginal) %>%
  ungroup() %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate_all(as.character)

log_debug("... split")
dataSplit <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferenceSplit,
    delim = "\t"
  ) %>%
  select(
    referenceOriginal = referenceOriginal_split,
    doi_split = referenceTranslatedDoi,
    journal_split = referenceTranslatedJournal,
    title_split = referenceTranslatedTitle,
    date_split = referenceTranslatedDate,
    author_split = referenceTranslatedAuthor,
    scoreCrossref_split = referenceTranslationScoreCrossref,
    scoreDistance_split = referenceTranslationScoreDistance
  ) %>%
  group_by(referenceOriginal) %>%
  mutate(level = row_number()) %>%
  relocate(level, .after = referenceOriginal) %>%
  ungroup() %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate_all(as.character)

log_debug("... full references")
dataFull <-
  read_delim(
    file = pathDataInterimTablesOriginalReferenceFull,
    delim = "\t",
    col_types = cols(.default = "c")
  )

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("...  cleaned organisms")
}
if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  dataCleanedOrganismManipulated_old <-
    read_delim(file = pathDataInterimDictionariesOrganismDictionary) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismType,
      organismValue,
      organismDetected
    ) %>%
    mutate_all(as.character)

  dataCleanedOrganismManipulated_new <-
    read_delim(file = pathDataInterimTablesCleanedOrganismFinal) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismType,
      organismValue,
      organismDetected
    ) %>%
    mutate_all(as.character)

  dataCleanedOrganismManipulated <-
    bind_rows(
      dataCleanedOrganismManipulated_old,
      dataCleanedOrganismManipulated_new
    ) %>%
    distinct()
}

if (!file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("... cleaned organisms")
}
if (!file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  dataCleanedOrganismManipulated <-
    read_delim(file = pathDataInterimTablesCleanedOrganismFinal) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismType,
      organismValue,
      organismDetected
    ) %>%
    mutate_all(as.character)
}

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  log_debug("... reference dictionary, this may take a while")
}
if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceDictionary <-
    read_delim(
      file = pathDataInterimDictionariesReferenceDictionary,
      delim = "\t",
      col_types = cols(.default = "c")
    )
}

log_debug("joining ...")
log_debug("... all reference types")
dataCrossref <- bind_rows(
  dataDoi,
  dataOriginal,
  dataPublishingDetails,
  dataPubmed,
  dataSplit,
  dataTitle
)

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  dataCrossref <- bind_rows(dataCrossref, referenceDictionary)
}

dataCrossref <- dataCrossref %>%
  filter(!is.na(referenceOriginal)) %>%
  distinct(
    referenceOriginal,
    referenceTranslatedType,
    origin,
    referenceTranslatedValue,
    level
  )

log_debug("... with organisms")
dataJoined <-
  left_join(dataFull, dataCleanedOrganismManipulated) %>%
  filter(!is.na(referenceValue)) %>%
  distinct(
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    organismDetected
  )

rm(
  dataDoi,
  dataOriginal,
  dataPublishingDetails,
  dataPubmed,
  dataSplit,
  dataTitle,
  dataFull,
  referenceDictionary,
  dataCleanedOrganismManipulated
)

log_debug("... with reference dictionary")
dataTranslated <- left_join(dataJoined,
  dataCrossref,
  by = c("referenceValue" = "referenceOriginal")
) %>%
  filter(
    !is.na(referenceTranslatedValue) |
      referenceType == "external" |
      referenceType == "journal" |
      referenceType == "authors" |
      referenceType == "isbn"
  ) %>%
  distinct(
    organismType,
    organismValue,
    organismDetected,
    referenceType,
    referenceValue,
    referenceTranslatedType,
    referenceTranslatedValue,
    level
  )

rm(dataJoined)

log_debug("ensuring directories exist")

ifelse(
  test = !dir.exists(pathDataInterimDictionaries),
  yes = dir.create(pathDataInterimDictionaries),
  no = paste(pathDataInterimDictionaries, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimDictionariesReference),
  yes = dir.create(pathDataInterimDictionariesReference),
  no = paste(pathDataInterimDictionariesReference, "exists")
)

log_debug("exporting, this may take a while if running full mode")
log_debug(pathDataInterimTablesTranslatedReferenceFile)
write_delim(
  x = dataTranslated,
  delim = "t",
  file = pathDataInterimTablesTranslatedReferenceFile
)

log_debug(pathDataInterimDictionariesReferenceDictionary)
write_delim(
  x = dataCrossref,
  delim = "t",
  file = pathDataInterimDictionariesReferenceDictionary
)

rm(dataCrossref)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
