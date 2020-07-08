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
reflistDoi <-
  pbmclapply(
    FUN = getrefDoi,
    X = dataDoi$referenceOriginalDoi,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE, 
    ignore.interactive = TRUE
  )

# joining with original dataframe
for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translatedDoi"] <-
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

for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translatedJournal"] <-
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

for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translatedTitle"] <-
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

for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translatedDate"] <-
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

for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translatedAuthor"] <-
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

for (i in 1:length(reflistDoi)) {
  dataDoi[i, "translationScore"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["score"]]),
        yes = reflistDoi[[i]][["data"]][["score"]],
        no = 0
      ),
      no = 0
    )[1])
}

dataDoi <- dataDoi %>% 
  mutate_all(as.character)
