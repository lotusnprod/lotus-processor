# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataSplit <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceSplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
if (nrow(dataSplit) != 1)
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit,
      X = dataSplit$referenceOriginal_split,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 2),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )

cat("This may take several minutes \n")

if (nrow(dataSplit) != 1)
  dataSplit <- getAllReferences(data = dataSplit,
                                referenceType = "split",
                                method = "osa")

if (nrow(dataSplit) == 1)
  dataSplit <- data.frame() %>%
  mutate(
    referenceOriginal_split = NA,
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
  x = dataSplit,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceSplit,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
