cat("This script integrates all reference translations together \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
source("r/vroom_safe.R")

cat("... files ... \n")
cat("... DOI \n")
dataDoi <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferenceDoi) %>%
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

cat("... original references \n")
dataOriginal <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferenceOriginal) %>%
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

cat("... PMID \n")
dataPubmed <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferencePubmed) %>%
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

cat("... publishing details \n")
dataPublishingDetails <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferencePublishingDetails) %>%
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

cat("... titles \n")
dataTitle <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferenceTitle) %>%
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

cat("... split \n")
dataSplit <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedReferenceSplit) %>%
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

cat("... full references \n")
dataFull <-
  vroom_read_safe(path = pathDataInterimTablesOriginalReferenceFull) %>%
  mutate_all(as.character)

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  cat("...  cleaned organisms \n")
}
if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  dataCleanedOrganismManipulated_old <-
    vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismOriginal,
      organismDetected
    ) %>%
    mutate_all(as.character)

  dataCleanedOrganismManipulated_new <-
    vroom_read_safe(path = pathDataInterimTablesCleanedOrganismFinal) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismOriginal,
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
  cat("... cleaned organisms \n")
}
if (!file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  dataCleanedOrganismManipulated <-
    vroom_read_safe(path = pathDataInterimTablesCleanedOrganismFinal) %>%
    mutate(
      organismDetected =
        word(organismDetected, 1)
    ) %>%
    distinct(
      organismOriginal,
      organismDetected
    ) %>%
    mutate_all(as.character)
}

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  cat("... reference dictionary, this may take a while \n")
}
if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesReferenceDictionary)
}

cat("joining ... \n")
cat("... all reference types \n")
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

cat("... with organisms \n")
dataJoined <-
  left_join(dataFull, dataCleanedOrganismManipulated) %>%
  filter(!is.na(referenceValue)) %>%
  distinct(
    organismOriginal,
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

cat("... with reference dictionary \n")
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
    organismOriginal,
    organismDetected,
    referenceType,
    referenceValue,
    referenceTranslatedType,
    referenceTranslatedValue,
    level
  )

rm(dataJoined)

cat("ensuring directories exist \n")

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

cat("exporting, this may take a while if running full mode \n")
cat(pathDataInterimTablesTranslatedReferenceFile, "\n")
vroom_write_safe(
  x = dataTranslated,
  path = pathDataInterimTablesTranslatedReferenceFile
)

cat(pathDataInterimDictionariesReferenceDictionary, "\n")
vroom_write_safe(
  x = dataCrossref,
  path = pathDataInterimDictionariesReferenceDictionary
)

rm(dataCrossref)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
