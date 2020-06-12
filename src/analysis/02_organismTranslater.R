# title: "Organism translatoR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

# loading files
## standard db
dataOrganism <- read_delim(
  file = gzfile(pathOriginalOrganism),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

## tcm names
tcmNamesDic <- read_delim(
  file = gzfile(pathInterimTcmNamesDic),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## common names
commonNamesDic <- read_delim(
  file = gzfile(pathInterimCommonNamesDic),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## taxa
problematicNamesDic <- read_delim(
  file = gzfile(pathInterimProblematicNamesDic),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

# manipulating organisms names
## removing some unnecessary characters
dataOrganism$organismTranslated <-
  trimws(dataOrganism$organismOriginal)

dataOrganism$organismTranslated <- gsub(
  pattern = "_",
  replacement = " ",
  x = dataOrganism$organismTranslated,
  fixed = TRUE
)

dataOrganism$organismTranslated <- gsub(
  pattern = " ",
  replacement = " ",
  x = dataOrganism$organismTranslated,
  fixed = TRUE
)

## replacing problematic names
a <- paste("\\b", problematicNamesDic$name, "\\b", sep = "")
b <- problematicNamesDic$nameMutated
X <- 1:nrow(dataOrganism)

replaceProblematicNames <- function(X) {
  dataOrganism$organismTranslated[X] <- stri_replace_all_regex(
    str = dataOrganism$organismTranslated[X],
    pattern = a,
    replacement = b,
    case_insensitive = FALSE,
    vectorize_all = FALSE
  )
}

dataOrganism$organismTranslated <- pbmclapply(
  FUN = replaceProblematicNames,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

## replacing tcm names
c <- paste("\\b", tcmNamesDic$vernacularName, "\\b", sep = "")
d <- tcmNamesDic$canonicalName

replaceTCMNames <- function(X) {
  dataOrganism$organismTranslated[X] <- stri_replace_all_regex(
    str = dataOrganism$organismTranslated[X],
    pattern = c,
    replacement = d,
    case_insensitive = TRUE,
    vectorize_all = FALSE
  )
}

dataOrganism$organismTranslated <- pbmclapply(
  FUN = replaceTCMNames,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

## replacing tcm names (2nd)
tcmNamesDic2 <- tcmNamesDic %>%
  filter(canonicalName != newCanonicalName)

e <-
  paste("\\b", tcmNamesDic2$canonicalName, "\\b", sep = "")
f <- tcmNamesDic2$newCanonicalName

replaceTCMNames2 <- function(X) {
  dataOrganism$organismTranslated[X] <- stri_replace_all_regex(
    str = dataOrganism$organismTranslated[X],
    pattern = e,
    replacement = f,
    case_insensitive = TRUE,
    vectorize_all = FALSE
  )
}

dataOrganism$organismTranslated <- pbmclapply(
  FUN = replaceTCMNames2,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

## replacing common names
g <- paste("\\b", commonNamesDic$vernacularName, "\\b", sep = "")
h <- commonNamesDic$canonicalName

replaceCommonNames <- function(X) {
  dataOrganism$organismTranslated[X] <- stri_replace_all_regex(
    str = dataOrganism$organismTranslated[X],
    pattern = g,
    replacement = h,
    case_insensitive = TRUE,
    vectorize_all = FALSE
  )
}

dataOrganism$organismTranslated <- pbmclapply(
  FUN = replaceCommonNames,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

## replacing problematic names (reverse)
i <- paste("\\b", problematicNamesDic$nameMutated, "\\b", sep = "")
j <- problematicNamesDic$name

replaceProblematicNamesReverse <- function(X) {
  dataOrganism$organismTranslated[X] <- stri_replace_all_regex(
    str = dataOrganism$organismTranslated[X],
    pattern = i,
    replacement = j,
    case_insensitive = FALSE,
    vectorize_all = FALSE
  )
}

dataOrganism$organismTranslated <- pbmclapply(
  FUN = replaceProblematicNamesReverse,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

dataOrganism$organismTranslated <-
  trimws(dataOrganism$organismTranslated)

# reorganizing
dataOrganism <- dataOrganism %>%
  select(organismOriginal,
         organismTranslated)

# in case of emergency for comparison
# dataOrganism <- dataOrganism %>%
#   select(organismOriginal) %>%
#   mutate(organismTranslated = organismOriginal)

# exporting
## bio
write.table(
  x = dataOrganism,
  file = gzfile(
    description = pathTranslatedOrganism,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
