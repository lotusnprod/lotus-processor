# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataTitle <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceTitle),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
if (nrow(dataTitle) != 1)
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit,
      X = dataTitle$referenceOriginal_title,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

if (nrow(dataTitle) != 1)
  dataTitle <- getBestReference(data = dataTitle,
                                referenceType = "title",
                                method = "osa")

if (nrow(dataTitle) == 1)
  dataTitle <- data.frame() %>%
  mutate(
    referenceOriginal_title = NA,
    referenceTranslatedDoi = NA,
    referenceTranslatedJournal = NA,
    referenceTranslatedTitle = NA,
    referenceTranslatedDate = NA,
    referenceTranslatedAuthor = NA,
    referenceTranslationScoreCrossref = NA,
    referenceTranslationScoreDistance = NA
  )

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

## exporting
write.table(
  x = dataTitle,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceTitle,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
