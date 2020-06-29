# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions/reference.R")

# joining all types together again
dataFull <- rbind(
  dataReferenceFillAuto,
  dataReferenceFillDoi,
  dataReferenceFillNoArticle,
  dataReferenceFillPubmed
)

# joining with original df
dataReferenceMin <- dataReferenceLongSplit %>%
  select(referenceOriginal)

dataReferenced <- left_join(dataReferenceMin, dataFull)

### problematic reference field containing multiple references with no clear association ###
### I also think we should output a "nonarticleref" column for external DB links etc ###
### inconsistency of journal name depending on retrieval method (check JNP) ###

dataReferencedSelected <- dataReferenced %>%
  select(
    referenceOriginal,
    referenceSplit = value,
    translatedDoi,
    translatedJournal,
    translatedTitle,
    translatedDate,
    translatedAuthor,
    translationScore
  ) %>%
  mutate(translationScore = replace_na(translationScore, "0")) %>%
  distinct(
    referenceOriginal,
    referenceSplit,
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

## ref
write.table(
  x = dataReferencedSelected,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReference,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
