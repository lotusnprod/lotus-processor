source("r/log_debug.R")
log_debug("This script performs original reference translation from crossRef")

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
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/progressr.R")

log_debug("loading original references lists")
length <-
  length(
    list.files(
      path = pathDataInterimTablesOriginalReferenceOriginalFolder,
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
create_dir_with_rm(export = pathDataInterimTablesTranslatedReferenceOriginalFolder)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferenceOriginalFolder,
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
      pathDataInterimTablesTranslatedReferenceOriginalFolder,
      "/",
      stringr::str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  dataOriginal <- readr::read_delim(
    file = inpath,
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataOriginal) != 0) {
    reflist <-
      getref_noLimit(xs = dataOriginal$referenceOriginal_original)
    log_debug("treating results, may take a while if full mode")
    dataOriginal2 <-
      getAllReferences(
        data = dataOriginal,
        referenceType = "original",
        method = "osa"
      )
  } else {
    dataOriginal2 <- data.frame() |>
      dplyr::mutate(
        referenceOriginal_original = NA,
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
    x = dataOriginal2,
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
dataOriginal3 <- do.call(
  "rbind",
  lapply(
    list.files(
      path = file.path(pathDataInterimTablesTranslatedReferenceOriginalFolder),
      pattern = "*.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      readr::read_delim(
        file = gzfile(
          description = file.path(pathDataInterimTablesTranslatedReferenceOriginalFolder, x)
        ),
        delim = "\t",
        col_types = cols(.default = "c"),
        trim_ws = TRUE,
        escape_double = TRUE,
      ) |>
        dplyr::mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferenceOriginal)
readr::write_delim(
  x = dataOriginal3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceOriginal,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
