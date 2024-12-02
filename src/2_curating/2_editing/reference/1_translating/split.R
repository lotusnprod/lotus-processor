source("r/log_debug.R")
log_debug("This script performs split references translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)
library(stringr)

log_debug("... functions")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")

log_debug("loading split references lists")
length <-
  length(
    list.files(
      path = pathDataInterimTablesOriginalReferenceSplitFolder,
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
create_dir_with_rm(export = pathDataInterimTablesTranslatedReferenceSplitFolder)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferenceSplitFolder,
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
      pathDataInterimTablesTranslatedReferenceSplitFolder,
      "/",
      stringr::str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  dataSplit <- readr::read_delim(
    file = inpath,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataSplit) != 0) {
    reflist <- getref_noLimit(xs = dataSplit$referenceOriginal_split) |>
      progressr::with_progress()
    log_debug("treating results, may take a while if full mode")
    dataSplit2 <-
      getAllReferences(
        data = dataSplit,
        referenceType = "split",
        method = "osa"
      )
  } else {
    dataSplit2 <- data.frame() |>
      dplyr::mutate(
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

  log_debug("exporting ...")
  readr::write_delim(
    x = dataSplit2,
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
dataSplit3 <- do.call(
  "rbind",
  lapply(
    list.files(
      path = file.path(pathDataInterimTablesTranslatedReferenceSplitFolder),
      pattern = "*.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      readr::read_delim(
        file = gzfile(
          description = file.path(pathDataInterimTablesTranslatedReferenceSplitFolder, x)
        ),
        delim = "\t",
        col_types = cols(.default = "c"),
        locale = locales,
        escape_double = TRUE,
        trim_ws = TRUE
      ) |>
        dplyr::mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferenceSplit)
readr::write_delim(
  x = dataSplit3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceSplit,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
