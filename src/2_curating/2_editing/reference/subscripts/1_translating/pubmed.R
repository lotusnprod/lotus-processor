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
reflistPubmed <- invisible(
  pbmclapply(
    FUN = getrefPubmed,
    X = as.character(dataPubmed$referenceOriginalPubmed),
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = 2,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

reflistPubmedBound <- bind_rows(reflistPubmed)

# joining with original dataframe
for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslatedDoi"] <-
    reflistPubmedBound[i, "translatedDoi"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslatedJournal"] <-
    reflistPubmedBound[i, "translatedJournal"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslatedTitle"] <-
    reflistPubmedBound[i, "translatedTitle"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslatedAuthor"] <-
    reflistPubmedBound[i, "translatedAuthor"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslatedDate"] <-
    reflistPubmedBound[i, "translatedDate"]
}

for (i in 1:nrow(reflistPubmedBound)) {
  dataPubmed[i, "referenceTranslationScore"] <- 1
}

dataPubmed <- dataPubmed %>%
  mutate_all(as.character)

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
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
