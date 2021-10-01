source("r/log_debug.R")
log_debug("This script performs publishing details translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(pbmcapply)
library(readr)
library(stringr)

log_debug("... functions")
source("r/getref_noLimit_publishingDetails.R")
source("r/getAllReferences.R")

log_debug("loading title lists")
length <-
  length(
    list.files(
      path = pathDataInterimTablesOriginalReferencePublishingDetailsFolder,
      pattern = "tsv"
    )
  )

cut <- 1000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

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

ifelse(
  test = !dir.exists(pathDataInterimTablesTranslatedReferencePublishingDetailsFolder),
  yes = dir.create(pathDataInterimTablesTranslatedReferencePublishingDetailsFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesTranslatedReferencePublishingDetailsFolder,
      full.names = TRUE
    )
  ) &
    dir.create(
      pathDataInterimTablesTranslatedReferencePublishingDetailsFolder,
      showWarnings = FALSE
    )
)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferencePublishingDetailsFolder,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv"
    )
  
  outpath <-
    paste0(
      pathDataInterimTablesTranslatedReferencePublishingDetailsFolder,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )
  
  log_debug(paste("step", i / cut, "of", length))
  
  dataPublishingDetails <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = TRUE,
    trim_ws = TRUE
  )
  
  log_debug("submitting to crossRef")
  if (nrow(dataPublishingDetails) != 0) {
    reflist <- invisible(
      pbmclapply(
        FUN = getref_noLimit_publishingDetails,
        X = dataPublishingDetails$referenceOriginal_publishingDetails,
        mc.preschedule = FALSE,
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
    log_debug("treating results, may take a while if full mode")
    dataPublishingDetails <-
      getAllReferences(
        data = dataPublishingDetails,
        referenceType = "publishingDetails",
        method = "osa"
      )
  } else {
    dataPublishingDetails2 <- data.frame() %>%
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
  
  log_debug("exporting ...")
  write_delim(
    x = dataPublishingDetails2,
    delim = "\t",
    file = outpath
  )
  
  ## cleaning memory
  gc(
    verbose = TRUE,
    reset = TRUE,
    full = TRUE
  )
}

log_debug("joining results with original lists")
dataPublishingDetails3 <- do.call(
  "rbind",
  lapply(
    list.files(
      path = file.path(pathDataInterimTablesTranslatedReferenceTitleFolder),
      pattern = "*.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      read_delim(
        file = gzfile(
          file.path(pathDataInterimTablesTranslatedReferencePublishingDetailsFolder, x)
        ),
        delim = "\t",
        escape_double = TRUE,
        trim_ws = TRUE
      ) %>%
        mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferencePublishingDetails)
write_delim(
  x = dataPublishingDetails3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferencePublishingDetails
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
