# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataPublishingDetails <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferencePublishingDetails),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
if (nrow(dataPublishingDetails) != 1)
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit,
      X = dataPublishingDetails$referenceOriginal_publishingDetails,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

print(x = "This may take several minutes")

# joining with original dataframe
if (nrow(dataPublishingDetails) != 1)
  dataPublishingDetails <-
  getAllReferences(data = dataPublishingDetails,
                   referenceType = "publishingDetails",
                   method = "osa")

if (nrow(dataPublishingDetails) == 1)
  dataPublishingDetails <- data.frame() %>%
  mutate(
    referenceOriginal_publishingDetails = NA,
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
  x = dataPublishingDetails,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferencePublishingDetails,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
