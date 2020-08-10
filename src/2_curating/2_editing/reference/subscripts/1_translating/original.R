# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

library(stringr)

## file
length <-
  length(
    list.files(path = pathDataInterimTablesOriginalReferenceOriginalFolder,
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
  !dir.exists(pathDataInterimTablesTranslatedReferenceOriginalFolder),
  dir.create(pathDataInterimTablesTranslatedReferenceOriginalFolder),
  FALSE
)

for (i in num) {
  inpath <- paste(
    pathDataInterimTablesOriginalReferenceOriginalFolder,
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
      pathDataInterimTablesTranslatedReferenceOriginalFolder,
      str_pad(
        string = i,
        width = 6  ,
        pad = "0"
      ),
      ".tsv.gz",
      sep = ""
    )
  
  print(paste("step", i / cut, "of", length))
  
  dataOriginal <- read_delim(
    file = inpath,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  # getting references
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
  
  print("treating results, may take a while if full mode")
  dataOriginal2 <-
    getAllReferences(data = dataOriginal,
                     referenceType = "original",
                     method = "osa")
  
  ## exporting
  write.table(
    x = dataOriginal2,
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

dataOriginal3 <- do.call("rbind",
                         lapply(list.files(
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
                             escape_double = FALSE,
                             trim_ws = TRUE
                           ) %>%
                             mutate_all(as.character)
                         }))
## exporting
write.table(
  x = dataOriginal3,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceOriginal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
