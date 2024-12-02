source("r/log_debug.R")
log_debug("This script performs DOI translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(readr)

log_debug("... functions")
source("r/getrefDoi.R")
source("r/progressr.R")

packageVersion("rcrossref")

log_debug("loading DOI list")
dataDoi <-
  readr::read_delim(
    file = pathDataInterimTablesOriginalReferenceDoi,
    delim = "\t",
    col_types = cols(.default = "c")
  )

log_debug("submitting to crossRef")
if (nrow(dataDoi) != 0) {
  reflistDoi <-
    getrefDoi(xs = dataDoi$referenceOriginal_doi)

  log_debug("joining results with original list")
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
} else {
  dataDoi <- data.frame() |>
    dplyr::mutate(
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
dataDoi <- dataDoi |>
  dplyr::mutate_all(as.character)

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
create_dir(export = pathDataInterimTablesTranslatedReference)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferenceDoi)
readr::write_delim(
  x = dataDoi,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceDoi,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
