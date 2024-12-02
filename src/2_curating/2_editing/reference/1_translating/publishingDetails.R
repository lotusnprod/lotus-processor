source("r/log_debug.R")
log_debug("This script performs publishing details translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(readr)
library(stringr)

log_debug("... functions")
source("r/getref_noLimit_publishingDetails.R")
source("r/getAllReferences.R")
source("r/progressr.R")

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
create_dir(export = pathDataInterimTablesTranslatedReference)
create_dir_with_rm(export = pathDataInterimTablesTranslatedReferencePublishingDetailsFolder)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferencePublishingDetailsFolder,
      "/",
      stringr::str_pad(
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
      stringr::str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  dataPublishingDetails <- readr::read_delim(
    file = inpath,
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataPublishingDetails) != 0) {
    reflist <- getref_noLimit_publishingDetails(xs = dataPublishingDetails$referenceOriginal_publishingDetails) |>
      progressr::with_progress()
    log_debug("treating results, may take a while if full mode")
    dataPublishingDetails2 <-
      getAllReferences(
        data = dataPublishingDetails,
        referenceType = "publishingDetails",
        method = "osa"
      )
  } else {
    dataPublishingDetails2 <- data.frame() |>
      dplyr::mutate(
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
  readr::write_delim(
    x = dataPublishingDetails2,
    delim = "\t",
    file = outpath,
    na = ""
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
      path = file.path(
        pathDataInterimTablesTranslatedReferencePublishingDetailsFolder
      ),
      pattern = "*.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      readr::read_delim(
        file = gzfile(
          description = file.path(
            pathDataInterimTablesTranslatedReferencePublishingDetailsFolder,
            x
          )
        ),
        delim = "\t",
        col_types = cols(.default = "c"),
        escape_double = TRUE,
        trim_ws = TRUE
      ) |>
        mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferencePublishingDetails)
readr::write_delim(
  x = dataPublishingDetails3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferencePublishingDetails,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
