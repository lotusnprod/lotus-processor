# title: "Ref integration"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## files
### doi
dataDoi <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceDoi),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

### original
dataOriginal <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

### pubmed
dataPubmed <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferencePubmed),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

### publishing details
dataPublishingDetails <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferencePublishingDetails),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

### title
dataTitle <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceTitle),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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
  mutate_all(as.character) %>%
  ungroup() %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = c("referenceTranslatedType", "origin"),
    names_sep = "_",
    values_to = "referenceTranslatedValue",
    values_drop_na = TRUE
  ) %>%
  mutate_all(as.character)

### split
dataSplit <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceSplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
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

### full
dataFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalReferenceFull),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

### cleaned
if (file.exists(pathDataInterimDictionariesOrganismDictionary))
  dataCleanedOrganismManipulated <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesOrganismDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  select(organismOriginal,
         organismCleaned) %>%
  mutate_all(as.character)

if (!file.exists(pathDataInterimDictionariesOrganismDictionary))
  dataCleanedOrganismManipulated <- read_delim(
    file = gzfile(description = pathDataInterimTablesCleanedOrganismFinal),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  select(organismOriginal,
         organismCleaned) %>%
  mutate_all(as.character)


cat("loading reference dictionary, this may take a while \n")

### dictionary
if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesReferenceDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

# joining all types together again
dataCrossref <- bind_rows(dataDoi,
                          dataOriginal,
                          dataPublishingDetails,
                          dataPubmed,
                          dataSplit,
                          dataTitle)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  dataCrossref <-  bind_rows(dataCrossref, referenceDictionary)

dataCrossref <- dataCrossref %>%
  filter(!is.na(referenceOriginal)) %>%
  distinct(
    referenceOriginal,
    referenceTranslatedType,
    origin,
    referenceTranslatedValue,
    level
  )

# joining full and cleaned organisms
dataJoined <-
  left_join(dataFull, dataCleanedOrganismManipulated) %>%
  filter(!is.na(referenceValue)) %>%
  distinct(organismOriginal,
           referenceType,
           referenceValue,
           organismCleaned)

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

dataTranslated <- left_join(dataJoined,
                            dataCrossref,
                            by = c("referenceValue" = "referenceOriginal")) %>%
  filter(
    !is.na(referenceTranslatedValue) |
      referenceType == "external" |
      referenceType == "journal" |
      referenceType == "authors" |
      referenceType == "isbn"
  ) %>%
  distinct(
    organismOriginal,
    organismCleaned,
    referenceType,
    referenceValue,
    referenceTranslatedType,
    referenceTranslatedValue,
    level
  )

## exporting
cat("exporting, this may take a while if running full mode \n")

write.table(
  x = dataTranslated,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceFile,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

### dictionary
ifelse(
  !dir.exists(pathDataInterimDictionariesReference),
  dir.create(pathDataInterimDictionariesReference),
  FALSE
)

write.table(
  x = dataCrossref,
  file = gzfile(
    description = pathDataInterimDictionariesReferenceDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
