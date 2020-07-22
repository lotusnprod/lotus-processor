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
reflist <- invisible(
  pbmclapply(
    FUN = getref,
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
for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslatedDoi"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["doi"]]),
        yes = reflist[[i]][["data"]][["doi"]],
        no = NA
      ),
      no = NA
    )[1])
}

for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslatedJournal"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["container.title"]]),
        yes = reflist[[i]][["data"]][["container.title"]],
        no = NA
      ),
      no = NA
    )[1])
}

for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslatedTitle"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["title"]]),
        yes = reflist[[i]][["data"]][["title"]],
        no = NA
      ),
      no = NA
    )[1])
}

for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslatedDate"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["issued"]]),
        yes = reflist[[i]][["data"]][["issued"]],
        no = NA
      ),
      no = NA
    )[1])
}

for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslatedAuthor"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["author"]][[1]][["family"]][1]),
        yes = reflist[[i]][["data"]][["author"]][[1]][["family"]][1],
        no = NA
      ),
      no = NA
    )[1])
  
}

for (i in 1:length(reflist)) {
  dataPublishingDetails[i, "referenceTranslationScore"] <-
    as.character(ifelse(
      test = !is.na(reflist[[i]]),
      yes = ifelse(
        test = !is.null(reflist[[i]][["data"]][["score"]]),
        yes = reflist[[i]][["data"]][["score"]],
        no = 0
      ),
      no = 0
    )[1])
}

dataPublishingDetails <- dataPublishingDetails %>%
  mutate_all(as.character)

dataPublishingDetails[] <-
  lapply(dataPublishingDetails, function(x)
    gsub("\r\n", " ", x))
dataPublishingDetails[] <-
  lapply(dataPublishingDetails, function(x)
    gsub("\r", " ", x))
dataPublishingDetails[] <-
  lapply(dataPublishingDetails, function(x)
    gsub("\n", " ", x))
dataPublishingDetails[] <-
  lapply(dataPublishingDetails, function(x)
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
