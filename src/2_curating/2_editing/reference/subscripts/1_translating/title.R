# title: "Ref translatoR"

# loading 
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataTitle <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceTitle),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
reflist <- invisible(
  pbmclapply(
    FUN = getref,
    X = dataTitle$referenceOriginalTitle,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE, 
    ignore.interactive = TRUE
  )
)

dataTitle <- dataTitle %>%
  mutate_all(as.character)

# joining with original dataframe
for (i in 1:length(reflist)) {
  dataTitle[i, "translatedDoi"] <-
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
  dataTitle[i, "translatedJournal"] <-
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
  dataTitle[i, "translatedTitle"] <-
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
  dataTitle[i, "translatedDate"] <-
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
  dataTitle[i, "translatedAuthor"] <-
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
  dataTitle[i, "translationScore"] <-
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
