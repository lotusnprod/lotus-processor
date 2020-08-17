# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
length <-
  length(
    list.files(path = pathDataInterimTablesOriginalReferenceTitleFolder,
               pattern = 'tsv')
  )

cut <- 1000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

# exporting
## creating directories if they do not exist
ifelse(
  !dir.exists(pathDataInterimTablesTranslated),
  dir.create(pathDataInterimTablesTranslated),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesTranslatedReference),
  dir.create(pathDataInterimTablesTranslatedReference),
  FALSE
)

ifelse(
  !dir.exists(pathDataInterimTablesTranslatedReferenceTitleFolder),
  dir.create(pathDataInterimTablesTranslatedReferenceTitleFolder),
  no = file.remove(
    list.files(path = pathDataInterimTablesTranslatedReferenceTitleFolder,
               full.names = TRUE)
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
      width = 6  ,
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
        width = 6  ,
        pad = "0"
      ),
      ".tsv.gz",
      sep = ""
    )
  
  print(paste("step", i / cut, "of", length))
  
  dataTitle <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  # getting references
  if (nrow(dataTitle) != 1)
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
  
  if (nrow(dataTitle) == 1)
    reflist <- list(NA)
  
  print("treating results, may take a while if full mode")
  if (nrow(dataTitle) != 0)
    dataTitle2 <-
    getAllReferences(data = dataTitle,
                     referenceType = "title",
                     method = "osa")
  
  if (nrow(dataTitle) == 0)
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
  
  ## exporting
  write.table(
    x = dataTitle2,
    file = gzfile(
      description = outpath,
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = TRUE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
  
  ## cleaning memory
  gc(verbose = TRUE,
     reset = TRUE,
     full = TRUE)
}

dataTitle3 <- do.call("rbind",
                      lapply(list.files(
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
                          escape_double = FALSE,
                          trim_ws = TRUE
                        ) %>%
                          mutate_all(as.character)
                      }))

## exporting
write.table(
  x = dataTitle3,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceTitle,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)