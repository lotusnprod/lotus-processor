# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataFull <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceFull),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character)

# joining all types together again
dataFull <- full_join(dataFull, dataDoi)
dataFull <- full_join(dataFull, dataPubmed)
dataFull <- full_join(dataFull, dataTitle)
dataFull <- full_join(dataFull, dataUnsplit)

### problematic reference field containing multiple references with no clear association ###
### I also think we should output a "nonarticleref" column for external DB links etc ###
### inconsistency of journal name depending on retrieval method (check JNP) ###

dataReferencedSelected <- dataFull %>%
  select(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    translatedDoi,
    translatedJournal,
    translatedTitle,
    translatedDate,
    translatedAuthor,
    translationScore
  ) %>%
  mutate(translationScore = replace_na(translationScore, "0")) %>%
  distinct(
    referenceOriginalAuthors,
    referenceOriginalDoi,
    referenceOriginalExternal,
    referenceOriginalIsbn,
    referenceOriginalJournal,
    referenceOriginalPubmed,
    referenceOriginalTitle,
    referenceOriginalUnsplit,
    translatedTitle,
    translatedJournal,
    translatedDate,
    translatedAuthor,
    translatedDoi,
    translationScore,
    .keep_all = TRUE
  )

dataReferencedSelected$translationScore <-
  as.numeric(dataReferencedSelected$translationScore)

dataReferencedSelected$translationScore[dataReferencedSelected$translationScore == 1] <-
  100

# exporting
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesTranslated),
  dir.create(pathDataInterimTablesTranslated),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesTranslatedReference),
  dir.create(pathDataInterimTablesTranslatedReference),
  FALSE
)

## ref
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
