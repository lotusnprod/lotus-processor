source("r/log_debug.R")
log_debug("This script performs original reference translation from crossRef")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(tidyverse)
library(pbmcapply)

log_debug("... functions")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/vroom_safe.R")

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
  test = !dir.exists(pathDataInterimTablesTranslatedReferenceOriginalFolder),
  yes = dir.create(pathDataInterimTablesTranslatedReferenceOriginalFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesTranslatedReferenceOriginalFolder,
      full.names = TRUE
    )
  ) &
    dir.create(
      pathDataInterimTablesTranslatedReferenceOriginalFolder,
      showWarnings = FALSE
    )
)

for (i in num) {
  inpath <- paste0(pathDataInterimTablesOriginalReferenceOriginalFolder, "/", str_pad(
    string = i,
    width = 6,
    pad = "0"
  ), ".tsv")

  outpath <-
    paste0(pathDataInterimTablesTranslatedReferenceOriginalFolder, "/", str_pad(
      string = i,
      width = 6,
      pad = "0"
    ), ".tsv.gz")

  log_debug(paste("step", i / cut, "of", length))

  dataOriginal <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataOriginal) != 1) {
    reflist <- invisible(
      pbmclapply(
        FUN = getref_noLimit,
        X = dataOriginal$referenceOriginal_original,
        mc.preschedule = FALSE,
        mc.set.seed = TRUE,
        mc.silent = TRUE,
        mc.cores = (parallel::detectCores() - 2),
        mc.cleanup = TRUE,
        mc.allow.recursive = TRUE,
        ignore.interactive = TRUE
      )
    )
  }

  log_debug("treating results, may take a while if full mode")
  if (nrow(dataOriginal) != 0) {
    dataOriginal2 <-
      getAllReferences(
        data = dataOriginal,
        referenceType = "original",
        method = "osa"
      )
  }

  if (nrow(dataOriginal) == 0) {
    dataOriginal2 <- data.frame() %>%
      mutate(
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
  vroom_write_safe(
    x = dataOriginal2,
    path = outpath
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
      read_delim(
        file = gzfile(
          file.path(pathDataInterimTablesTranslatedReferenceOriginalFolder, x)
        ),
        delim = "\t",
        trim_ws = TRUE,
        escape_double = TRUE,
      ) %>%
        mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesTranslatedReferenceOriginal)
vroom_write_safe(
  x = dataOriginal3,
  path = pathDataInterimTablesTranslatedReferenceOriginal
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
