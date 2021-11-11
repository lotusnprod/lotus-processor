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
  test = !dir.exists(pathDataInterimTablesTranslatedReferenceSplitFolder),
  yes = dir.create(pathDataInterimTablesTranslatedReferenceSplitFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesTranslatedReferenceSplitFolder,
      full.names = TRUE
    )
  ) &
    dir.create(
      pathDataInterimTablesTranslatedReferenceSplitFolder,
      showWarnings = FALSE
    )
)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesOriginalReferenceSplitFolder,
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
      pathDataInterimTablesTranslatedReferenceSplitFolder,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  dataSplit <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = TRUE,
    trim_ws = TRUE
  )

  log_debug("submitting to crossRef")
  if (nrow(dataSplit) != 0) {
    reflist <- invisible(
      pbmclapply(
        FUN = getref_noLimit,
        X = dataSplit$referenceOriginal_split,
        mc.preschedule = TRUE,
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
    dataSplit2 <-
      getAllReferences(
        data = dataSplit,
        referenceType = "split",
        method = "osa"
      )
  } else {
    dataSplit2 <- data.frame() %>%
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

  log_debug("exporting ...")
  write_delim(
    x = dataSplit2,
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
dataSplit3 <- do.call(
  "rbind",
  lapply(
    list.files(
      path = file.path(pathDataInterimTablesTranslatedReferenceSplitFolder),
      pattern = "*.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      read_delim(
        file = gzfile(
          file.path(pathDataInterimTablesTranslatedReferenceSplitFolder, x)
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
log_debug(pathDataInterimTablesTranslatedReferenceSplit)
write_delim(
  x = dataSplit3,
  delim = "\t",
  file = pathDataInterimTablesTranslatedReferenceSplit
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
