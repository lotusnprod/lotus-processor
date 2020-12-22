cat("This script performs taxonomy cleaning and alignment. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/y_as_na.R")
source("r/manipulating_taxo.R")
source("r/taxo_cleaning_auto.R")
source("r/rank_order.R")

cat("loading ... \n")
cat("... libraries \n")
library(tidyverse)

log_debug(" Step 4")
cat("... files ... \n")
cat("... cleaned organisms \n")
dataCleanedOrganism <- read_delim(
  file = gzfile(pathDataInterimTablesCleanedOrganismTranslatedTable),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = FALSE
) %>%
  distinct()

cat(" ... taxa ranks dictionary \n")
taxaRanksDictionary <- read_delim(
  file = pathDataInterimDictionariesTaxaRanks,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("manipulating taxonomic levels \n")
if (nrow(dataCleanedOrganism) != 0) {
  dataCleanedOrganismManipulated <-
    manipulating_taxo(
      dfsel = dataCleanedOrganism,
      dic = taxaRanksDictionary
    )
}

if (nrow(dataCleanedOrganism) == 0) {
  dataCleanedOrganismManipulated <- data.frame() %>%
    mutate(
      organismOriginal = NA,
      organismDetected = NA,
      organismCleaned = NA,
      organismCleanedId = NA,
      organismCleanedRank = NA,
      organismDbTaxo = NA,
      organismDbTaxoQuality = NA,
      organismTaxonIds = NA,
      organismTaxonRanks = NA,
      organismTaxonomy = NA,
      organism_1_kingdom = NA,
      organism_2_phylum = NA,
      organism_3_class = NA,
      organism_4_order = NA,
      organism_5_family = NA,
      organism_6_genus = NA,
      organism_7_species = NA,
      organism_8_variety = NA,
      # organism_1_kingdom_id = NA,
      # organism_2_phylum_id = NA,
      # organism_3_class_id = NA,
      # organism_4_order_id = NA,
      # organism_5_family_id = NA,
      # organism_6_genus_id = NA,
      # organism_7_species_id = NA,
      # organism_8_variety_id = NA
    )
}

dataCleanedOrganismManipulated <-
  dataCleanedOrganismManipulated[order(match(
    dataCleanedOrganismManipulated$organismCleanedRank,
    rank_order
  )), ]

dataCleanedOrganismManipulated_clean <- dataCleanedOrganismManipulated %>%
  distinct(organismOriginal, organismDetected, organismCleaned, .keep_all = TRUE) %>%
  select(organismOriginal, organismDetected, organismCleaned, organismCleanedRank)

dataCleanedOrganismManipulated_clean_2 <- left_join(
  dataCleanedOrganismManipulated_clean,
  dataCleanedOrganismManipulated
) %>% distinct()

dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCleanedOrganismManipulated_clean_2)

cat("selecting \n")
dataCuratedOrganismAuto[setdiff(
  x = c(
    "organismOriginal",
    "organismDetected",
    "organismCleaned",
    "organismCleanedId",
    "organismCleanedRank",
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
    "organism_8_variety"
    # "organism_1_kingdom_id",
    # "organism_2_phylum_id",
    # "organism_3_class_id",
    # "organism_4_order_id",
    # "organism_5_family_id",
    # "organism_6_genus_id",
    # "organism_7_species_id",
    # "organism_8_variety_id"
  ),
  y = names(dataCuratedOrganismAuto)
)] <- NA

dataCuratedOrganismAuto <- dataCuratedOrganismAuto %>%
  select(
    organismOriginal,
    organismDetected,
    organismCleaned,
    organismCleanedId,
    organismCleanedRank,
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
    # organism_1_kingdom_id,
    # organism_2_phylum_id,
    # organism_3_class_id,
    # organism_4_order_id,
    # organism_5_family_id,
    # organism_6_genus_id,
    # organism_7_species_id,
    # organism_8_variety_id,
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
