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
reflist <- invisible(
  pbmclapply(
    FUN = getref,
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

print(x = "This may take several minutes")

# joining with original dataframe
for (i in 1:length(reflist)) {
  dataSplit[i, "referenceTranslatedDoi"] <-
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
  dataSplit[i, "referenceTranslatedJournal"] <-
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
  dataSplit[i, "referenceTranslatedTitle"] <-
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
  dataSplit[i, "referenceTranslatedDate"] <-
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
  dataSplit[i, "referenceTranslatedAuthor"] <-
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
  dataSplit[i, "referenceTranslationScore"] <-
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

dataSplit <- dataSplit %>%
  mutate_all(as.character)

dataSplit[] <-
  lapply(dataSplit, function(x)
    gsub("\r\n", " ", x))
dataSplit[] <-
  lapply(dataSplit, function(x)
    gsub("\r", " ", x))
dataSplit[] <-
  lapply(dataSplit, function(x)
    gsub("\n", " ", x))
dataSplit[] <-
  lapply(dataSplit, function(x)
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
