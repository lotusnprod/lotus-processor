# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataDoi <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceDoi),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
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

# joining with original dataframe
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
