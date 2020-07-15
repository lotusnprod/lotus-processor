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


dataPubmed <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferencePubmed),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataPubmed)[2:ncol(dataPubmed)] <-
  paste(colnames(dataPubmed)[2:ncol(dataPubmed)], "pubmed", sep = "_")

dataTitle <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceTitle),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataTitle)[2:ncol(dataTitle)] <-
  paste(colnames(dataTitle)[2:ncol(dataTitle)], "title", sep = "_")

dataUnsplit <- read_delim(
  file = gzfile(pathDataInterimTablesTranslatedReferenceUnsplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

colnames(dataUnsplit)[2:ncol(dataUnsplit)] <-
  paste(colnames(dataUnsplit)[2:ncol(dataUnsplit)], "unsplit", sep = "_")

dataFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  distinct(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit
  ) %>%
  mutate_all(as.character)

# joining all types together again
dataFullWide <- full_join(dataFull, dataDoi)
dataFullWide <- full_join(dataFullWide, dataPubmed)
dataFullWide <- full_join(dataFullWide, dataTitle)
dataFullWide <- full_join(dataFullWide, dataUnsplit)

dataFullLong <- dataFullWide %>%
  pivot_longer(
    cols = 9:ncol(.),
    names_to = c(".value", "level"),
    names_sep = "_",
    values_to = "reference",
    values_drop_na = TRUE
  )

dataFullLongFilled <- left_join(dataFull, dataFullLong)

### inconsistency of journal name depending on retrieval method (check JNP) ###

dataReferencedSelected <- dataFullLongFilled %>%
  select(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    referenceTranslatedDoi,
    referenceTranslatedJournal,
    referenceTranslatedTitle,
    referenceTranslatedDate,
    referenceTranslatedAuthor,
    referenceTranslationScore
  ) %>%
  mutate(referenceTranslationScore = replace_na(referenceTranslationScore, "0")) %>%
  distinct(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    referenceTranslatedTitle,
    referenceTranslatedJournal,
    referenceTranslatedDate,
    referenceTranslatedAuthor,
    referenceTranslatedDoi,
    referenceTranslationScore,
    .keep_all = TRUE
  )

dataReferencedSelected$referenceTranslationScore <-
  as.numeric(dataReferencedSelected$referenceTranslationScore)

dataReferencedSelected$referenceTranslationScore[dataReferencedSelected$referenceTranslationScore == 1] <-
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
