# title: "Ref translatoR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# getting references
reflistDoi <- invisible(
  pbmclapply(
    FUN = getrefDoi,
    X = dataReferenceFillDoi$value,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

# joining with original dataframe
for (i in 1:length(reflistDoi)) {
  dataReferenceFillDoi[i, "translatedDoi"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["doi"]]),
        yes = reflistDoi[[i]][["data"]][["doi"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedJournal"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["container.title"]]),
        yes = reflistDoi[[i]][["data"]][["container.title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedTitle"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["title"]]),
        yes = reflistDoi[[i]][["data"]][["title"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedDate"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["issued"]]),
        yes = reflistDoi[[i]][["data"]][["issued"]],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translatedAuthor"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1]),
        yes = reflistDoi[[i]][["data"]][["author"]][[1]][["family"]][1],
        no = NA
      ),
      no = NA
    )[1])
  
  dataReferenceFillDoi[i, "translationScore"] <-
    as.character(ifelse(
      test = !is.na(reflistDoi[[i]][["data"]][["score"]]),
      yes = ifelse(
        test = !is.null(reflistDoi[[i]][["data"]][["score"]]),
        yes = reflistDoi[[i]][["data"]][["score"]],
        no = 0
      ),
      no = 0
    )[1])
}
