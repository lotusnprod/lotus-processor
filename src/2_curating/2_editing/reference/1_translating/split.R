source("r/log_debug.R")
log_debug("This script performs split references translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(readr)

log_debug("... functions")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")

log_debug("loading split references list")
dataSplit <-
  read_delim(
    file = pathDataInterimTablesOriginalReferenceSplit,
    delim = "\t"
  )

log_debug("submitting to crossRef")
if (!is.na(dataSplit[, 1])) {
  reflist <- invisible(
    pbmclapply(
      FUN = getref_noLimit,
      X = dataSplit$referenceOriginal_split,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE,
      mc.silent = TRUE,
      mc.cores = min(max(1, parallel::detectCores() - 1), 10),
      mc.cleanup = TRUE,
      mc.allow.recursive = TRUE,
      ignore.interactive = TRUE,
      mc.style = "txt",
      mc.substyle = 1
    )
  )

  log_debug("This may take several minutes")
  dataSplit <- getAllReferences(
    data = dataSplit,
    referenceType = "split",
    method = "osa"
  )
} else {
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
log_debug(pathDataInterimTablesTranslatedReferenceSplit)
write_delim(
  x = dataSplit,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceSplit
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
