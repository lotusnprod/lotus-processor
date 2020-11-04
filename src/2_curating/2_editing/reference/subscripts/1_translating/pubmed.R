cat("This script performs PMID translation from pubmed API \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/reference.R")

cat("loading PMID list \n")
dataPubmed <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferencePubmed),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references ##getting them with pubmed API and not crossRef because crossRef pubmed ID not working!!
## 2
# mc cores set to 2 because fails otherwise (entrez limitation probably)
cat("submitting to entrez \n")
if (nrow(dataPubmed) != 1) {
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
}
if (nrow(dataPubmed) != 1) {
  reflistPubmedBound <- bind_rows(reflistPubmed)
}

cat("joining results with original list \n")
if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedDoi"] <-
      reflistPubmedBound[i, "translatedDoi"]
  }
}

if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedJournal"] <-
      reflistPubmedBound[i, "translatedJournal"]
  }
}

if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedTitle"] <-
      reflistPubmedBound[i, "translatedTitle"]
  }
}

if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedAuthor"] <-
      reflistPubmedBound[i, "translatedAuthor"]
  }
}

if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslatedDate"] <-
      reflistPubmedBound[i, "translatedDate"]
  }
}

if (nrow(dataPubmed) != 1) {
  for (i in 1:nrow(reflistPubmedBound)) {
    dataPubmed[i, "referenceTranslationScoreCrossref"] <- 1
    dataPubmed[i, "referenceTranslationScoreDistance"] <- 0
  }
}

if (nrow(dataPubmed) == 1) {
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
}

cat("removing unfriendly characters \n")
dataPubmed <- dataPubmed %>%
  mutate_all(as.character)

dataPubmed[] <-
  lapply(dataPubmed, function(x) {
    gsub("\r\n", " ", x)
  })
dataPubmed[] <-
  lapply(dataPubmed, function(x) {
    gsub("\r", " ", x)
  })
dataPubmed[] <-
  lapply(dataPubmed, function(x) {
    gsub("\n", " ", x)
  })
dataPubmed[] <-
  lapply(dataPubmed, function(x) {
    gsub("\t", " ", x)
  })

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesTranslated),
  yes = dir.create(pathDataInterimTablesTranslated),
  no = paste(pathDataInterimTablesTranslated, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesTranslatedReference),
  yes = dir.create(pathDataInterimTablesTranslatedReference),
  no = paste(pathDataInterimTablesTranslatedReference, "exists")
)

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedReferencePubmed, "\n")
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")