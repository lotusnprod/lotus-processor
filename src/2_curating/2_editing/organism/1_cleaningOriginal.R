cat("This script performs canonical name recognition on the original organism field. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("...paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/gnfinder_cleaning.R")
source("r/vroom_safe.R")
source("r/split_data_table.R")
source("r/y_as_na.R")

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(tidyverse)

log_debug("  Step 1")
cat("... taxa ranks dictionary \n")
taxaRanksDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesTaxaRanks)

wrongVerifiedDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesTaxaWrongVerified) %>%
  as.list()

organismTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalOrganismFull) %>%
  distinct()

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

cat("submitting to GNVerifier \n")
system(command = paste("bash", pathOriginalGnverifierScript))

verified <-
  stream_in(con = file(pathDataInterimTablesCleanedOrganismVerifiedOriginalTable))

verified_df <- verified %>%
  data.frame() %>%
  select(
    -curation,
    -matchType
  ) %>%
  unnest(preferredResults, names_repair = "minimal") %>%
  filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") %>%
  filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) %>%
  mutate(organismType = "clean") %>%
  select(
    organismType,
    organismValue = input,
    organismCleaned = currentCanonicalFull,
    organismDbTaxo = dataSourceTitleShort,
    taxonId = currentRecordId,
    organismCleanedCurrent = currentName,
    taxonomy = classificationPath,
    rank = classificationRanks
  )

## example ID 165 empty, maybe fill later on
verified_df$organismDbTaxo <-
  y_as_na(verified_df$organismDbTaxo, "")

dataOrganismVerified <- left_join(
  organismTable,
  verified_df
) %>%
  select(
    organismType,
    organismValue,
    organismCleaned,
    organismDbTaxo,
    taxonId,
    organismCleanedCurrent,
    taxonomy,
    rank
  )

dataOrganismNoVerified <- dataOrganismVerified %>%
  arrange(organismDbTaxo) %>%
  distinct(organismValue, .keep_all = TRUE) %>%
  filter(is.na(organismDbTaxo)) %>%
  distinct(organismValue) %>%
  data.table()

dataOrganismVerified <- dataOrganismVerified %>%
  filter(!is.na(organismDbTaxo))

cat(pathDataInterimTablesOriginalOrganism, "\n")

if (nrow(dataOrganismNoVerified) != 0) {
  split_data_table(
    x = dataOrganismNoVerified,
    no_rows_per_frame = 2000,
    # else verification takes too long and gnverifier times out
    text = "",
    path_to_store = pathDataInterimTablesOriginalOrganism
  )
}

cat("submitting to GNFinder \n")
system(command = paste("bash", pathOriginalGnfinderScript))

cat("treating GNFinder results \n")
length <-
  length(list.files(
    path = pathDataInterimTablesOriginalOrganism,
    pattern = "^[0-9]{6}.tsv"
  ))

cut <-
  2000 # else verification takes too long and gnverifier times out

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
    cat(paste("step", j, "of", length, "\n"))
    tryCatch(
      {
        dataCleanOriginalOrganism[[j]] <-
          gnfinder_cleaning(
            num = i,
            organismCol = "organismValue"
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
      organismValue,
      organismCleaned = canonicalname,
      organismCleanedCurrent = canonicalnameCurrent,
      organismDbTaxo = dbTaxo,
      everything()
    ) %>%
    select(-ids, -dbQuality)
}

if (length(dataCleanOriginalOrganism) == 0) {
  dataCleanedOriginalOrganism <- data.frame() %>%
    mutate(
      organismValue = NA,
      organismCleaned = NA,
      organismCleanedCurrent = NA,
      organismDbTaxo = NA,
      taxonId = NA,
      taxonomy = NA,
      rank = NA,
    )
}

dataCleanedOriginalOrganismUnique <- dataCleanedOriginalOrganism %>%
  distinct(organismValue, organismCleaned, .keep_all = TRUE)

cat("exporting ... \n")
if (length != 0) {
  cat(pathDataInterimTablesCleanedOrganismOriginalTable, "\n")
}

if (length != 0) {
  vroom_write_safe(
    x = dataCleanedOriginalOrganism,
    path = pathDataInterimTablesCleanedOrganismOriginalTable
  )
}

if (length != 0) {
  cat(
    pathDataInterimTablesCleanedOrganismOriginalUniqueTable,
    "\n"
  )
}

if (length != 0) {
  vroom_write(
    x = dataCleanedOriginalOrganismUnique,
    path = gzfile(
      description = pathDataInterimTablesCleanedOrganismOriginalUniqueTable,
      compression = 9,
      encoding = "UTF-8"
    ),
    num_threads = 1,
    bom = TRUE,
    quote = "none",
    escape = "double",
    delim = "\t",
    col_names = TRUE,
    progress = TRUE,
    append = FALSE
  )
  ## because of univocity parser settings
}

vroom_write_safe(
  x = dataOrganismVerified,
  path = pathDataInterimTablesCleanedOrganismOriginalVerifiedTable
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
