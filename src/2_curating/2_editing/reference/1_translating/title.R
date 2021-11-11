source("r/log_debug.R")
log_debug("This script performs title translation from crossRef")

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
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/parallel.R")

log_debug("loading title lists")
length <-
  length(
    list.files(
      path = pathDataInterimTablesOriginalReferenceTitleFolder,
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
  test = !dir.exists(pathDataInterimTablesTranslatedReferenceTitleFolder),
  yes = dir.create(pathDataInterimTablesTranslatedReferenceTitleFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesTranslatedReferenceTitleFolder,
      full.names = TRUE
    )
  ) &
    dir.create(
      pathDataInterimTablesTranslatedReferenceTitleFolder,
      showWarnings = FALSE
    )
)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferenceTitleFolder,
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
      pathDataInterimTablesTranslatedReferenceTitleFolder,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  dataTitle <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataTitle) != 0) {
    reflist <- invisible(
      pbmclapply(
        FUN = getref_noLimit,
        X = dataTitle$referenceOriginal_title,
        mc.preschedule = FALSE,
        mc.set.seed = TRUE,
        mc.silent = TRUE,
        mc.cores = numCores,
        mc.cleanup = TRUE,
        mc.allow.recursive = TRUE,
        ignore.interactive = TRUE,
        mc.style = "txt",
        mc.substyle = 1
      )
    )
    log_debug("treating results, may take a while if full mode")
    dataTitle2 <-
      getAllReferences(
        data = dataTitle,
        referenceType = "title",
        method = "osa"
      )
  } else {
    dataTitle2 <- data.frame() %>%
      mutate(
        referenceOriginal_title = NA,
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
    x = dataTitle2,
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
dataTitle3 <- do.call(
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
          file.path(pathDataInterimTablesTranslatedReferenceTitleFolder, x)
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
log_debug(pathDataInterimTablesTranslatedReferenceTitle)
write_delim(
  x = dataTitle3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceTitle
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
