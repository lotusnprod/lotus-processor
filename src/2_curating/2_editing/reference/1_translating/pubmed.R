source("r/log_debug.R")
log_debug("This script performs PMID translation from pubmed API")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(readr)

log_debug("... functions")
source("r/getrefPubmed.R")

log_debug("loading PMID list")
dataPubmed <-
  read_delim(
    file = pathDataInterimTablesOriginalReferencePubmed,
    delim = "\t"
  )

# getting references ##getting them with pubmed API and not crossRef because crossRef pubmed ID not working!!
# mc cores set to 1 because fails otherwise (entrez limitation of 10 calls per sec probably)
log_debug("submitting to entrez")
if (!is.na(dataPubmed[, 1])) {
  reflistPubmed <- invisible(
    pbmclapply(
      FUN = getrefPubmed,
      X = as.character(dataPubmed$referenceOriginal_pubmed),
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = 1,
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
  )

  if (is.null(reflistPubmed$value)) {
    reflistPubmedBound <- bind_rows(reflistPubmed)
  } else {
    reflistPubmedBound <- bind_rows(reflistPubmed$value)
  }

  log_debug("joining results with original list")
  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslatedDoi"] <-
      reflistPubmedBound[i, "translatedDoi"]
  }

  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslatedJournal"] <-
      reflistPubmedBound[i, "translatedJournal"]
  }

  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslatedTitle"] <-
      reflistPubmedBound[i, "translatedTitle"]
  }

  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslatedAuthor"] <-
      reflistPubmedBound[i, "translatedAuthor"]
  }

  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslatedDate"] <-
      reflistPubmedBound[i, "translatedDate"]
  }

  for (i in seq_len(nrow(reflistPubmedBound))) {
    dataPubmed[i, "referenceTranslationScoreCrossref"] <- 1
    dataPubmed[i, "referenceTranslationScoreDistance"] <- 0
  }
} else {
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

log_debug("removing unfriendly characters")
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

log_debug("ensuring directories exist")
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

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferencePubmed)
write_delim(
  x = dataPubmed,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferencePubmed
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
