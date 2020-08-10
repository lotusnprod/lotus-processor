# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataSplit <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceSplit),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# getting references
reflist <- invisible(
  pbmclapply(
    FUN = getref_noLimit,
    X = dataSplit$referenceOriginal_split,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

print(x = "This may take several minutes")
dataSplit <- getAllReferences(data = dataSplit,
                              referenceType = "split",
                              method = "osa")

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

## exporting
write.table(
  x = dataSplit,
  file = gzfile(
    description = pathDataInterimTablesTranslatedReferenceSplit,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
