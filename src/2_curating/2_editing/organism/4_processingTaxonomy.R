source("r/log_debug.R")
log_debug("This script performs taxonomy cleaning and alignment. ")

start <- Sys.time()

log_debug("sourcing ... ")
log_debug("... paths ")
source("paths.R")

log_debug("... functions ")
source("r/y_as_na.R")
source("r/manipulating_taxo.R")
source("r/taxo_cleaning_auto.R")
source("r/rank_order.R")

log_debug("loading ... ")
log_debug("... libraries ")
library(dplyr)
library(jsonlite)
library(readr)
library(tidyr)

log_debug(" Step 4")
log_debug("... files ... ")
log_debug("... cleaned organisms ")
dataCleanedOrganism <-
  read_delim(
    file = pathDataInterimTablesProcessedOrganismTranslatedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  distinct()

log_debug(" ... taxa ranks dictionary ")
taxaRanksDictionary <-
  read_delim(
    file = pathDataInterimDictionariesTaxaRanks,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )
wrongVerifiedDictionary <-
  read_delim(
    file = pathDataInterimDictionariesTaxaWrongVerified,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  as.list()

dataCleanedOrganismVerify <- dataCleanedOrganism %>%
  filter(!is.na(organismCleaned)) %>%
  distinct(organismCleaned)

write_delim(
  x = dataCleanedOrganismVerify,
  file = gzfile(
    description = pathDataInterimTablesProcessedOrganismVerifyTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  na = "",
  delim = "\t",
  quote = "none",
  escape = "double"
)
## because gnverifier does not parse quotes

log_debug("submitting to GNVerifier")
if (.Platform$OS.type == "unix") {
  system(command = paste("bash", pathGnverifierScript))
} else {
  shell(paste("bash", pathGnverifierScript))
}

verified <-
  stream_in(con = file(pathDataInterimTablesProcessedOrganismVerifiedTable))

verified_df <- verified %>%
  data.frame() %>%
  select(-curation, -matchType) %>%
  unnest(results, names_repair = "minimal") %>%
  filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") %>%
  filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) %>%
  arrange(desc(sortScore)) %>%
  distinct(name, dataSourceTitleShort, .keep_all = TRUE) %>%
  select(
    organismCleaned = name,
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

log_debug("manipulating taxonomic levels ")
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

log_debug("selecting ")
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

log_debug("exporting ... ")
log_debug(pathDataInterimTablesProcessedOrganismFinal)
write_delim(
  x = dataCuratedOrganismAuto,
  delim = "\t",
  file = pathDataInterimTablesProcessedOrganismFinal,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
