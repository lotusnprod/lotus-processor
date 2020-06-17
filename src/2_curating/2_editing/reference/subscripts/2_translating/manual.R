# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# temporary
dataReferenceFillManual <- dataReferenceFillManual %>%
  mutate(
    translatedDoi = NA,
    translatedJournal = NA,
    translatedTitle = NA,
    translatedDate = NA,
    translatedAuthor = NA,
    translationScore = 0
  )
