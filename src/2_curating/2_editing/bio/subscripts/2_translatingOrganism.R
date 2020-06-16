# title: "treating bio"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

## tcm names
tcmNamesDic <- read_delim(
  file = gzfile(pathDataInterimDictionariesTcmNames),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## common names
commonNamesDic <- read_delim(
  file = gzfile(pathDataInterimDictionariesCommonNames),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
)

## black listed strings
blacklistDictionary <- read_delim(
  file = pathDataInterimDictionariesCommonBlackDic,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate(n = str_count(string = blackName)) %>%
  arrange(desc(n)) %>%
  select(-n)

dataInterimOrganism <- dataCleanedOriginalOrganism %>%
  mutate(
    organismInterim = stri_replace_all_regex(
      str =  organismOriginal,
      pattern = paste("\\b", blacklistDictionary$blackName, "\\b", sep = ""),
      replacement = "",
      case_insensitive = TRUE,
      vectorize_all = FALSE
    )
  ) %>%
  rowwise() %>%
  mutate(organismInterim = stri_replace_all_regex(
    str =  organismInterim,
    pattern = paste("\\b", organismCleaned, "\\b", sep = ""),
    replacement = "",
  ))

dataInterimOrganism$organismInterim <-
  gsub(
    pattern = ".",
    replacement = "",
    x = dataInterimOrganism$organismInterim,
    fixed = TRUE
  )

dataInterimOrganism$organismInterim <-
  y_as_na(x = dataInterimOrganism$organismInterim, y = "")

dataInterimOrganism$organismInterim <-
  trimws(dataInterimOrganism$organismInterim)

dataInterimOrganismToFill <- dataInterimOrganism %>%
  filter(!is.na(organismInterim))

# manipulating organisms names
## removing some unnecessary characters
dataInterimOrganismToFill$organismInterim <- gsub(
  pattern = "_",
  replacement = " ",
  x = dataInterimOrganismToFill$organismInterim,
  fixed = TRUE
)

dataInterimOrganismToFill$organismInterim <- gsub(
  pattern = " ",
  replacement = " ",
  x = dataInterimOrganismToFill$organismInterim,
  fixed = TRUE
)

## replacing tcm names
a <- paste("\\b", tcmNamesDic$vernacularName, "\\b", sep = "")
b <- tcmNamesDic$canonicalName
X <- 1:nrow(dataInterimOrganismToFill)

replaceTCMNames <- function(X) {
  dataInterimOrganismToFill$organismInterim[X] <-
    stri_replace_all_regex(
      str = dataInterimOrganismToFill$organismInterim[X],
      pattern = a,
      replacement = b,
      case_insensitive = TRUE,
      vectorize_all = FALSE
    )
}

dataInterimOrganismToFill$organismInterim <- pbmclapply(
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

c <-
  paste("\\b", tcmNamesDic2$canonicalName, "\\b", sep = "")
d <- tcmNamesDic2$newCanonicalName

replaceTCMNames2 <- function(X) {
  dataInterimOrganismToFill$organismInterim[X] <-
    stri_replace_all_regex(
      str = dataInterimOrganismToFill$organismInterim[X],
      pattern = c,
      replacement = d,
      case_insensitive = TRUE,
      vectorize_all = FALSE
    )
}

dataInterimOrganismToFill$organismInterim <- pbmclapply(
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
e <- paste("\\b", commonNamesDic$vernacularName, "\\b", sep = "")
f <- commonNamesDic$canonicalName

replaceCommonNames <- function(X) {
  dataInterimOrganismToFill$organismInterim[X] <-
    stri_replace_all_regex(
      str = dataInterimOrganismToFill$organismInterim[X],
      pattern = e,
      replacement = f,
      case_insensitive = TRUE,
      vectorize_all = FALSE
    )
}

dataInterimOrganismToFill$organismInterim <- pbmclapply(
  FUN = replaceCommonNames,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

dataInterimOrganismToFill$organismInterim <-
  as.character(dataInterimOrganismToFill$organismInterim)

dataInterimOrganismToFill$organismInterim <-
  y_as_na(x = dataInterimOrganismToFill$organismInterim, y = "")

dataInterimOrganismToFill$organismInterim <-
  y_as_na(x = dataInterimOrganismToFill$organismInterim, y = "NA")

dataInterimOrganismToFillGnfinder <- dataInterimOrganismToFill %>%
  data.table() %>%
  mutate_all(as.character) %>%
  filter(!is.na(organismInterim)) %>%
  select(organismInterim)

# exporting
split_data_table(
  x = dataInterimOrganismToFillGnfinder,
  no_rows_per_frame = cut,
  text = "translatedOrganismGnfinderUntil_",
  path_to_store = pathTranslatedOrganismDistinct
)
