source("r/log_debug.R")
log_debug("This script performs DOI translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(readr)

log_debug("... functions")
source("r/getrefDoi.R")

log_debug("loading DOI list")
dataDoi <-
  read_delim(
    file = pathDataInterimTablesOriginalReferenceDoi,
    delim = "\t"
  )

log_debug("submitting to crossRef")
if (nrow(dataDoi) != 1) {
  reflistDoi <-
    pbmclapply(
      FUN = getrefDoi,
      X = dataDoi$referenceOriginal_doi,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 1),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
}

log_debug("joining results with original list")
if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslatedDoi"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["doi"]]),
          yes = reflistDoi[[i]][["data"]][["doi"]],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslatedJournal"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["container.title"]]),
          yes = reflistDoi[[i]][["data"]][["container.title"]],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslatedTitle"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["title"]]),
          yes = reflistDoi[[i]][["data"]][["title"]],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslatedDate"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["issued"]]),
          yes = reflistDoi[[i]][["data"]][["issued"]],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslatedAuthor"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1]),
          yes = reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslationScoreCrossref"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["score"]]),
          yes = reflistDoi[[i]][["data"]][["score"]],
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) != 1) {
  for (i in seq_along(reflistDoi)) {
    dataDoi[i, "referenceTranslationScoreDistance"] <-
      as.character(ifelse(
        test = !is.na(reflistDoi[[i]]),
        yes = ifelse(
          test = !is.null(reflistDoi[[i]][["data"]][["doi"]]),
          yes = 0,
          no = NA
        ),
        no = NA
      )[1])
  }
}

if (nrow(dataDoi) == 1) {
  dataDoi <- data.frame() %>%
    mutate(
      referenceOriginal_doi = NA,
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
dataDoi <- dataDoi %>%
  mutate_all(as.character)

dataDoi[] <-
  lapply(dataDoi, function(x) {
    gsub("\r\n", " ", x)
  })
dataDoi[] <-
  lapply(dataDoi, function(x) {
    gsub("\r", " ", x)
  })
dataDoi[] <-
  lapply(dataDoi, function(x) {
    gsub("\n", " ", x)
  })
dataDoi[] <-
  lapply(dataDoi, function(x) {
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
log_debug(pathDataInterimTablesTranslatedReferenceDoi)
write_delim(
  x = dataDoi,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceDoi
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
