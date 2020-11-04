cat("This script performs canonical name recognition on the original organism field. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("...paths \n")
source("paths.R")

cat("... functions \n")
source("functions/bio.R")
source("functions/helpers.R")
source("functions/log.R")
source("2_curating/2_editing/organism/functions/manipulating_taxo.R") # shouldnt these path be in the path.R ???
source("2_curating/2_editing/organism/functions/gnfinder_cleaning.R") # shouldnt these path be in the path.R ???

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(jsonlite)
library(tidyverse)

log_debug("  Step 1")
cat("... taxa ranks dictionary \n")
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesCleaned),
  yes = dir.create(pathDataInterimTablesCleaned),
  no = paste(pathDataInterimTablesCleaned, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesCleanedOrganism),
  yes = dir.create(pathDataInterimTablesCleanedOrganism),
  no = paste(pathDataInterimTablesCleanedOrganism, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesCleanedOrganismOriginal),
  yes = dir.create(pathDataInterimTablesCleanedOrganismOriginal),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesCleanedOrganismOriginal,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesCleanedOrganismOriginal,
      showWarnings = FALSE
    )
)

cat("submitting to GNFinder \n")
system(command = paste("bash", pathOriginalGnfinderScript))

length <-
  length(list.files(
    path = pathDataInterimTablesOriginalOrganism,
    pattern = "tsv"
  ))

cut <- 10000

if (length != 0) {
  num <- as.integer(seq(
    from = 1 * cut,
    to = length * cut,
    by = cut
  ))
}

dataCleanOriginalOrganism <- list()

cat("cleaning GNFinder output \n")
if (length != 0) {
  for (i in num) {
    j <- i / cut
    cat(paste("step", j, "of", length))
    tryCatch(
      {
        dataCleanOriginalOrganism[[j]] <-
          gnfinder_cleaning(
            num = i,
            organismCol = "organismOriginal"
          )
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }
}

cat("selecting and reordering \n")
if (length(dataCleanOriginalOrganism) != 0) {
  dataCleanedOriginalOrganism <-
    bind_rows(dataCleanOriginalOrganism) %>%
    select(
      organismOriginal,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      everything()
    ) %>%
    select(-nchar, -sum)
}

if (length(dataCleanOriginalOrganism) == 0) {
  dataCleanedOriginalOrganism <- data.frame() %>%
    mutate(
      organismOriginal = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxo = NA,
      value_min = NA,
      value_max = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
      ids = NA,
      dbQuality = NA
    )
}

dataCleanedOriginalOrganismUnique <-
  dataCleanedOriginalOrganism %>%
  distinct(organismOriginal, organismCleaned, .keep_all = TRUE)

cat("exporting ... \n")
if (length != 0) {
  cat(pathDataInterimTablesCleanedOrganismOriginalTable, "\n")
}

if (length != 0) {
  write.table(
    x = dataCleanedOriginalOrganism,
    file = gzfile(
      description = pathDataInterimTablesCleanedOrganismOriginalTable,
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (length != 0) {
  cat(
    pathDataInterimTablesCleanedOrganismOriginalUniqueTable,
    "\n"
  )
}

if (length != 0) {
  write.table(
    x = dataCleanedOriginalOrganismUnique,
    file = gzfile(
      description = pathDataInterimTablesCleanedOrganismOriginalUniqueTable,
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")