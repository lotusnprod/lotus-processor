cat("This script performs DOI translation from crossRef \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/reference.R")

cat("loading DOI list \n")
dataDoi <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceDoi),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("submitting to crossRef \n")
if (nrow(dataDoi) != 1)
  reflistDoi <-
  pbmclapply(
    FUN = getrefDoi,
    X = dataDoi$referenceOriginal_doi,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )

cat("joining results with original list \n")
if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) != 1)
  for (i in 1:length(reflistDoi)) {
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

if (nrow(dataDoi) == 1)
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

cat("removing unfriendly characters \n")
dataDoi <- dataDoi %>%
  mutate_all(as.character)

dataDoi[] <-
  lapply(dataDoi, function(x)
    gsub("\r\n", " ", x))
dataDoi[] <-
  lapply(dataDoi, function(x)
    gsub("\r", " ", x))
dataDoi[] <-
  lapply(dataDoi, function(x)
    gsub("\n", " ", x))
dataDoi[] <-
  lapply(dataDoi, function(x)
    gsub("\t", " ", x))

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
cat(pathDataInterimTablesTranslatedReferenceDoi, "\n")
write.table(
  x = dataDoi,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceDoi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

end <- Sys.time()

cat("Script finished in", end - start , "seconds \n")
