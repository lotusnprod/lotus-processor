cat("This script performs split references translation from crossRef \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(pbmcapply)

cat("... functions \n")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/vroom_safe.R")

cat("loading split references list \n")
dataSplit <-
  vroom_read_safe(path = pathDataInterimTablesOriginalReferenceSplit)

cat("submitting to crossRef \n")
if (nrow(dataSplit) != 1) {
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit,
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
}

cat("This may take several minutes \n")

if (nrow(dataSplit) != 1) {
  dataSplit <- getAllReferences(
    data = dataSplit,
    referenceType = "split",
    method = "osa"
  )
}

if (nrow(dataSplit) == 1) {
  dataSplit <- data.frame() %>%
    mutate(
      referenceOriginal_split = NA,
      referenceTranslatedDoi = NA,
      referenceTranslatedJournal = NA,
      referenceTranslatedTitle = NA,
      referenceTranslatedDate = NA,
      referenceTranslatedAuthor = NA,
      referenceTranslationScoreCrossref = NA,
      referenceTranslationScoreDistance = NA
    )
}

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
cat(pathDataInterimTablesTranslatedReferenceSplit, "\n")
vroom_write_safe(
  x = dataSplit,
  path = pathDataInterimTablesTranslatedReferenceSplit
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
