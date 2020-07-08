# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataUnsplit <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceUnsplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
reflist <- invisible(
  pbmclapply(
    FUN = getref,
    X = dataUnsplit$referenceOriginalUnsplit,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

# joining with original dataframe
for (i in 1:length(reflist)) {
  dataUnsplit[i, "translatedDoi"] <-
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
  dataUnsplit[i, "translatedJournal"] <-
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
  dataUnsplit[i, "translatedTitle"] <-
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
  dataUnsplit[i, "translatedDate"] <-
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
  dataUnsplit[i, "translatedAuthor"] <-
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
  dataUnsplit[i, "translationScore"] <-
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

dataUnsplit <- dataUnsplit %>%
  mutate_all(as.character)
