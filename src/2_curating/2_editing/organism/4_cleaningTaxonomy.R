cat("This script performs taxonomy cleaning and alignment. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/y_as_na.R")
source("r/taxo_cleaning_auto.R")

cat("loading ... \n")
cat("... libraries \n")
library(tidyverse)

log_debug(" Step 4")
cat("... files ... \n")
cat("... cleaned organisms \n")
dataCleanedOrganismManipulated <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  select(everything(), organismDetected = organismCleaned) %>%
  relocate(organismDetected, .after = organismTaxonomy) %>%
  distinct()

dataCuratedOrganism <- dataCleanedOrganismManipulated %>%
  mutate(organismCleaned = as.character(apply(dataCleanedOrganismManipulated[7:15], 1, function(x) {
    tail(na.omit(x), 1)
  })))

dataCuratedOrganism$organismCleaned <-
  y_as_na(
    x = dataCuratedOrganism$organismCleaned,
    y = "character(0)"
  )

dataCuratedOrganism$organismCleaned <-
  y_as_na(
    x = dataCuratedOrganism$organismCleaned,
    y = "NA"
  )

dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCuratedOrganism)

cat("selecting \n")
dataCuratedOrganismAuto[setdiff(
  x = c(
    "organismOriginal",
    "organismCleaned",
    "organismDbTaxo",
    "organismDbTaxoQuality",
    "organismTaxonIds",
    "organismTaxonRanks",
    "organismTaxonomy",
    "organism_1_kingdom",
    "organism_2_phylum",
    "organism_3_class",
    "organism_4_order",
    "organism_5_family",
    "organism_6_genus",
    "organism_7_species",
    "organism_8_quality"
  ),
  y = names(dataCuratedOrganismAuto)
)] <- NA

dataCuratedOrganismAuto <- dataCuratedOrganismAuto %>%
  select(
    organismOriginal,
    organismDetected,
    organismCleaned,
    organismDbTaxo,
    organismDbTaxoQuality,
    organismTaxonIds,
    organismTaxonRanks,
    organismTaxonomy,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    organism_7_species,
    organism_8_variety
  ) %>%
  filter(grepl(pattern = "[[:alnum:]]", x = organismTaxonRanks))

cat("exporting ... \n")
cat(pathDataInterimTablesCleanedOrganismFinal, "\n")
write.table(
  x = dataCuratedOrganismAuto,
  file = gzfile(
    description = pathDataInterimTablesCleanedOrganismFinal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
