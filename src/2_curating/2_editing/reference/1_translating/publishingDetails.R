source("r/log_debug.R")
log_debug("This script performs publishing details translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(tidyr)

log_debug("... functions")
source("r/getref_noLimit_publishingDetails.R")
source("r/getAllReferences.R")

log_debug("loading publishing details list")
dataPublishingDetails <-
  read_delim(
    file = pathDataInterimTablesOriginalReferencePublishingDetails,
    delim = "\t"
  )

log_debug("submitting to crossRef")
if (nrow(dataPublishingDetails) != 1) {
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit_publishingDetails,
      X = dataPublishingDetails$referenceOriginal_publishingDetails,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = (parallel::detectCores() - 1),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
  )
}

log_debug("This may take several minutes")

log_debug("joining results with original list")
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

log_debug("ensuring directories exist")
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

log_debug("exporting ...")
log_debug(
  pathDataInterimTablesTranslatedReferencePublishingDetails
)
write_delim(
  x = dataPublishingDetails,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferencePublishingDetails
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
