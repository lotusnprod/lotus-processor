cat("This script performs original reference translation from crossRef \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(pbmcapply)

cat("... functions \n")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/vroom_safe.R")

cat("loading original references lists \n")
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
  inpath <- paste(
    pathDataInterimTablesOriginalReferenceOriginalFolder,
    str_pad(
      string = i,
      width = 6,
      pad = "0"
    ),
    ".tsv",
    sep = ""
  )

  outpath <-
    paste(
      pathDataInterimTablesTranslatedReferenceOriginalFolder,
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz",
      sep = ""
    )

  cat(paste("step", i / cut, "of", length))

  dataOriginal <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

  cat("submitting to crossRef \n")
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

  cat("treating results, may take a while if full mode \n")
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

  cat("exporting ... \n")
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

cat("joining results with original lists \n")
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
        quote = ""
      ) %>%
        mutate_all(as.character)
    }
  )
)

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedReferenceOriginal, "\n")
vroom_write_safe(
  x = dataOriginal3,
  path = pathDataInterimTablesTranslatedReferenceOriginal
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")