# title: "Ref translatoR"

# loading
## paths
source("paths.R")

## functions
source("functions/reference.R")

## file
dataOriginal <- read_delim(
  file = gzfile(pathDataInterimTablesOriginalReferenceOriginal),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# splitting because of too big results otherwise
## 1
dataOriginal_1 <- dataOriginal %>%
  filter(row_number() %% 2 == 0)

## 2
dataOriginal_2 <- dataOriginal %>%
  filter(row_number() %% 2 == 1)

dataOriginal <- rbind(dataOriginal_1, dataOriginal_2)

# getting references
## 1
reflist_1 <- invisible(
  pbmclapply(
    FUN = getref_noLimit,
    X = dataOriginal_1$referenceOriginal_original,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

## 2
reflist_2 <- invisible(
  pbmclapply(
    FUN = getref_noLimit,
    X = dataOriginal_2$referenceOriginal_original,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

## joining results
reflist <- append(reflist_1, reflist_2)

print(x = "This may take several minutes")

# joining with original dataframe
dataOriginal <-
  getAllReferences(data = dataOriginal,
                   referenceType = "original",
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
  x = dataOriginal,
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
