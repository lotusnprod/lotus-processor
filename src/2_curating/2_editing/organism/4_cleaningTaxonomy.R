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
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... libraries \n")
library(tidyverse)
library(jsonlite)

log_debug(" Step 4")
cat("... files ... \n")
cat("... cleaned organisms \n")
dataCleanedOrganism <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismTranslatedTable) %>%
  distinct()

cat(" ... taxa ranks dictionary \n")
taxaRanksDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesTaxaRanks)

dataCleanedOrganismVerify <- dataCleanedOrganism %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(organismCleaned)

vroom_write(
  x = dataCleanedOrganismVerify,
  path = gzfile(
    description = pathDataInterimTablesCleanedOrganismVerifyTable,
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
## because gnverify does not parse quotes

cat("submitting to GNVerify \n")
system(command = paste("bash", pathGnverifyScript))

verified <-
  stream_in(con = file(pathDataInterimTablesCleanedOrganismVerifiedTable))

verified_df <- verified %>%
  data.frame() %>%
  select(-curation, -matchType) %>%
  unnest(preferredResults, names_repair = "minimal") %>%
  filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") %>%
  select(
    organismCleaned = input,
    organismDbTaxo = dataSourceTitleShort,
    taxonId = currentRecordId,
    currentName,
    currentCanonicalFull,
    taxonomy = classificationPath,
    rank = classificationRanks
  )

## example ID 165 empty, maybe fill later on
verified_df$organismDbTaxo <-
  y_as_na(verified_df$organismDbTaxo, "")

dataCleanedOrganismVerified <- left_join(
  dataCleanedOrganism,
  verified_df
) %>%
  select(
    organismType,
    organismValue,
    organismCleaned,
    organismDbTaxo,
    taxonId,
    currentName,
    currentCanonicalFull,
    taxonomy,
    rank
  ) %>%
  filter(!is.na(organismDbTaxo))

## last version does not contain "species"
indexFungorum <- dataCleanedOrganismVerified %>%
  filter(organismDbTaxo == "Index Fungorum") %>%
  mutate(
    rank = ifelse(
      test = rank == "kingdom|phylum|class|order|family|",
      yes = "kingdom|phylum|class|order|family|species",
      no = rank
    )
  )

if (nrow(indexFungorum != 0)) {
  dataCleanedOrganismVerified <- dataCleanedOrganismVerified %>%
    filter(organismDbTaxo != "Index Fungorum") %>%
    bind_rows(., indexFungorum)
}

cat("manipulating taxonomic levels \n")
if (nrow(dataCleanedOrganism) != 0) {
  dataCleanedOrganismManipulated <-
    manipulating_taxo(
      dfsel = dataCleanedOrganismVerified,
      dic = taxaRanksDictionary
    )
}

if (nrow(dataCleanedOrganism) == 0) {
  dataCleanedOrganismManipulated <- data.frame() %>%
    mutate(
      organismType = NA,
      organismValue = NA,
      organismDetected = NA,
      organismCleaned = NA,
      organismCleanedId = NA,
      organismCleanedRank = NA,
      organismDbTaxo = NA,
      organismDbTaxoQuality = NA,
      name = NA,
      organismTaxonIds = NA,
      organismTaxonRanks = NA,
      organismTaxonomy = NA,
      organism_1_kingdom = NA,
      organism_2_phylum = NA,
      organism_3_class = NA,
      organism_4_order = NA,
      organism_5_family = NA,
      organism_6_genus = NA,
      organism_6_1_subgenus = NA,
      organism_7_species = NA,
      organism_7_1_subspecies = NA,
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

dataCleanedOrganismManipulated_clean <-
  dataCleanedOrganismManipulated %>%
  distinct(
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    .keep_all = TRUE
  ) %>%
  select(
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    organismCleanedRank
  )

dataCleanedOrganismManipulated_clean_2 <-
  left_join(
    dataCleanedOrganismManipulated_clean,
    dataCleanedOrganismManipulated
  ) %>%
  distinct()

dataCuratedOrganismAuto <-
  taxo_cleaning_auto(dfsel = dataCleanedOrganismManipulated_clean_2)

cat("selecting \n")
dataCuratedOrganismAuto[setdiff(
  x = c(
    "organismType",
    "organismValue",
    "organismDetected",
    "organismCleaned",
    "organismCleanedId",
    "organismCleanedRank",
    "organismDbTaxo",
    "organismDbTaxoQuality",
    # "organismTaxonIds",
    "organismTaxonRanks",
    "organismTaxonomy",
    "organism_1_kingdom",
    "organism_2_phylum",
    "organism_3_class",
    "organism_4_order",
    "organism_5_family",
    "organism_6_genus",
    "organism_6_1_subgenus",
    "organism_7_species",
    "organism_7_1_subspecies",
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
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    organismCleanedId,
    organismCleanedRank,
    organismDbTaxo,
    # organismDbTaxoQuality,
    # organismTaxonIds,
    organismTaxonRanks,
    organismTaxonomy,
    organism_1_kingdom,
    organism_2_phylum,
    organism_3_class,
    organism_4_order,
    organism_5_family,
    organism_6_genus,
    # organism_6_1_subgenus,
    organism_7_species,
    # organism_7_1_subspecies,
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
vroom_write_safe(
  x = dataCuratedOrganismAuto,
  path = pathDataInterimTablesCleanedOrganismFinal
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
