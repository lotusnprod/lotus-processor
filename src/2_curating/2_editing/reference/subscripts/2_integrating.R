# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## files
dataDoi <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceDoi),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataDoi)[2:ncol(dataDoi)] <-
  paste(colnames(dataDoi)[2:ncol(dataDoi)], "doi", sep = "_")

dataOriginal <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataOriginal)[2:ncol(dataOriginal)] <-
  paste(colnames(dataOriginal)[2:ncol(dataOriginal)], "original", sep = "_")

dataPubmed <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferencePubmed),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataPubmed)[2:ncol(dataPubmed)] <-
  paste(colnames(dataPubmed)[2:ncol(dataPubmed)], "pubmed", sep = "_")

dataPublishingDetails <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferencePublishingDetails),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataPublishingDetails)[2:ncol(dataPublishingDetails)] <-
  paste(colnames(dataPublishingDetails)[2:ncol(dataPublishingDetails)], "publishingDetails", sep = "_")

dataTitle <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceTitle),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataTitle)[2:ncol(dataTitle)] <-
  paste(colnames(dataTitle)[2:ncol(dataTitle)], "title", sep = "_")

dataSplit <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceSplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataSplit)[2:ncol(dataSplit)] <-
  paste(colnames(dataSplit)[2:ncol(dataSplit)], "split", sep = "_")

dataFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  distinct(
    organismOriginal,
    referenceOriginal_authors,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_isbn,
    referenceOriginal_journal,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_publishingDetails,
    referenceOriginal_title,
    referenceOriginal_split
  ) %>%
  mutate_all(as.character)

# joining all types together again
dataFullWide <- full_join(dataFull, dataDoi)
dataFullWide <- full_join(dataFullWide, dataPubmed)
dataFullWide <- full_join(dataFullWide, dataTitle)
dataFullWide <- full_join(dataFullWide, dataPublishingDetails)
dataFullWide <- full_join(dataFullWide, dataSplit)
dataFullWide <- full_join(dataFullWide, dataOriginal)


dataFullLong <- dataFullWide %>%
  pivot_longer(
    cols = 12:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "reference",
    values_drop_na = TRUE
  )

dataFullLongFilled <- left_join(dataFull, dataFullLong)

dataReferencedSelected <- dataFullLongFilled %>%
  select(
    organismOriginal,
    referenceOriginal_authors,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_isbn,
    referenceOriginal_journal,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_publishingDetails,
    referenceOriginal_title,
    referenceOriginal_split,
    referenceTranslatedDoi,
    referenceTranslatedJournal,
    referenceTranslatedTitle,
    referenceTranslatedDate,
    referenceTranslatedAuthor,
    referenceTranslationScoreCrossref,
    referenceTranslationScoreDistance
  ) %>%
  distinct(
    organismOriginal,
    referenceOriginal_authors,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_isbn,
    referenceOriginal_journal,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_publishingDetails,
    referenceOriginal_title,
    referenceOriginal_split,
    referenceTranslatedTitle,
    referenceTranslatedJournal,
    referenceTranslatedDate,
    referenceTranslatedAuthor,
    referenceTranslatedDoi,
    referenceTranslationScoreCrossref,
    referenceTranslationScoreDistance,
    .keep_all = TRUE
  )

dataReferencedSelected$referenceTranslationScoreCrossref <-
  as.numeric(dataReferencedSelected$referenceTranslationScoreCrossref)

dataReferencedSelected$referenceTranslationScoreDistance <-
  as.numeric(dataReferencedSelected$referenceTranslationScoreDistance)

dataReferencedSelected$referenceTranslationScoreCrossref[dataReferencedSelected$referenceTranslationScoreCrossref == 1] <-
  100

## exporting
write.table(
  x = dataReferencedSelected,
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
