# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataPubmed <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferencePubmed),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references ##getting them with pubmed API and not crossRef because crossRef pubmed ID not working!!
## 2
# mc cores set to 2 because fails otherwise (entrez limitation probably)
if (nrow(dataPubmed) != 1)
  reflistPubmed <- invisible(
    pbmclapply(
      FUN = getrefPubmed,
      X = as.character(dataPubmed$referenceOriginal_pubmed),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = 2,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE
    )
  )
if (nrow(dataPubmed) != 1)
  reflistPubmedBound <- bind_rows(reflistPubmed)

# joining with original dataframe
if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedDoi"] <-
      reflistPubmedBound[i, "translatedDoi"]
  }

if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedJournal"] <-
      reflistPubmedBound[i, "translatedJournal"]
  }

if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedTitle"] <-
      reflistPubmedBound[i, "translatedTitle"]
  }

if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedAuthor"] <-
      reflistPubmedBound[i, "translatedAuthor"]
  }

if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedDate"] <-
      reflistPubmedBound[i, "translatedDate"]
  }

if (nrow(dataPubmed) != 1)
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslationScoreCrossref"] <- 1
    dataPubmed[i, "referenceTranslationScoreDistance"] <- 0
  }

if (nrow(dataPubmed) == 1)
  dataPubmed <- data.frame() %>%
  mutate(
    referenceOriginal_pubmed = NA,
    referenceTranslatedDoi = NA,
    referenceTranslatedJournal = NA,
    referenceTranslatedTitle = NA,
    referenceTranslatedDate = NA,
    referenceTranslatedAuthor = NA,
    referenceTranslationScoreCrossref = NA,
    referenceTranslationScoreDistance = NA
  )

dataPubmed <- dataPubmed %>%
  mutate_all(as.character)

dataPubmed[] <-
  lapply(dataPubmed, function(x)
    gsub("\r\n", " ", x))
dataPubmed[] <-
  lapply(dataPubmed, function(x)
    gsub("\r", " ", x))
dataPubmed[] <-
  lapply(dataPubmed, function(x)
    gsub("\n", " ", x))
dataPubmed[] <-
  lapply(dataPubmed, function(x)
    gsub("\t", " ", x))

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
  x = dataPubmed,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferencePubmed,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
