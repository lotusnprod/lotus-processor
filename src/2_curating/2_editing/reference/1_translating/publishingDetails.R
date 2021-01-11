cat("This script performs publishing details translation from crossRef \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(pbmcapply)

cat("... functions \n")
source("r/getref_noLimit_publishingDetails.R")
source("r/getAllReferences.R")
source("r/vroom_safe.R")

cat("loading publishing details list \n")
dataPublishingDetails <-
  vroom_read_safe(path = pathDataInterimTablesOriginalReferencePublishingDetails)

cat("submitting to crossRef \n")
if (nrow(dataPublishingDetails) != 1) {
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit_publishingDetails,
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
}

cat("This may take several minutes \n")

cat("joining results with original list \n")
if (nrow(dataPublishingDetails) != 1) {
  dataPublishingDetails <-
    getAllReferences(
      data = dataPublishingDetails,
      referenceType = "publishingDetails",
      method = "osa"
    )
}

if (nrow(dataPublishingDetails) == 1) {
  dataPublishingDetails <- data.frame() %>%
    mutate(
      referenceOriginal_publishingDetails = NA,
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
cat(
  pathDataInterimTablesTranslatedReferencePublishingDetails,
  "\n"
)
vroom_write_safe(
  x = dataPublishingDetails,
  path = pathDataInterimTablesTranslatedReferencePublishingDetails
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
