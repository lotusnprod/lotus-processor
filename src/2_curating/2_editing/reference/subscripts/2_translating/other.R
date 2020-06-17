# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")


# selecting fields corresponding to other sources than articles
dataReferenceFillNoArticle <- dataReferenceLongSplit %>%
  filter(
    grepl(pattern = "^KNApSAcK Database$",
          x = value) |
      grepl(pattern = "^DrDuke$",
            x = value) |
      grepl(pattern = "^patent:",
            x = value) |
      grepl(pattern = "^CHEBI:",
            x = value)
  ) %>%
  mutate(
    translatedDoi = NA,
    translatedJournal = NA,
    translatedTitle = NA,
    translatedDate = NA,
    translatedAuthor = NA,
    translationScore = 0
  )
