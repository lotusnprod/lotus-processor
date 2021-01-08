cat("This script performs title translation from crossRef \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
groundhog.library(tidyverse, date = groundhog.day)
groundhog.library(pbmcapply, date = groundhog.day)

cat("... functions \n")
source("r/getref_noLimit.R")
source("r/getAllReferences.R")
source("r/vroom_safe.R")

cat("loading title lists \n")
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
  inpath <- paste(
    pathDataInterimTablesOriginalReferenceTitleFolder,
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
      pathDataInterimTablesTranslatedReferenceTitleFolder,
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv.gz",
      sep = ""
    )

  cat(paste("step", i / cut, "of", length))

  dataTitle <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

  cat("submitting to crossRef \n")
  if (nrow(dataTitle) != 1) {
    reflist <- invisible(
      pbmclapply(
        FUN = getref_noLimit,
        X = dataTitle$referenceOriginal_title,
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

  if (nrow(dataTitle) == 1) {
    reflist <- list(NA)
  }

  cat("treating results, may take a while if full mode \n")
  if (nrow(dataTitle) != 0) {
    dataTitle2 <-
      getAllReferences(
        data = dataTitle,
        referenceType = "title",
        method = "osa"
      )
  }

  if (nrow(dataTitle) == 0) {
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

  cat("exporting ... \n")
  vroom_write_safe(
    x = dataTitle2,
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
        quote = "",
        trim_ws = TRUE
      ) %>%
        mutate_all(as.character)
    }
  )
)

cat("exporting ... \n")
cat(pathDataInterimTablesTranslatedReferenceTitle, "\n")
vroom_write_safe(
  x = dataTitle3,
  path = pathDataInterimTablesTranslatedReferenceTitle
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
