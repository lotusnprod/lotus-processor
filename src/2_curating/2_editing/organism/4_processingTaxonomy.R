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
log_debug("... files ...")
log_debug("... cleaned organisms")
dataCleanedOrganism <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismTranslatedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct()

log_debug("... taxa ranks dictionary")
taxaRanksDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaRanks,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )
log_debug("... manual fixes dictionaries")
wrongVerifiedDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaWrongVerified,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  as.list()

wrongHomonymsDictionary <-
  readr::read_delim(
    file = pathDataInterimDictionariesTaxaWrongHomonyms,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

dataCleanedOrganismVerify <- dataCleanedOrganism |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::distinct(organismCleaned)

readr::write_delim(
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
  jsonlite::stream_in(con = file(pathDataInterimTablesProcessedOrganismVerifiedTable))

verified_df <- verified |>
  data.frame() |>
  dplyr::select(-curation, -matchType) |>
  tidyr::unnest(results, names_repair = "minimal") |>
  dplyr::filter(dataSourceTitleShort != "IRMNG (old)" &
    dataSourceTitleShort != "IPNI") |>
  dplyr::filter(!matchedName %in% wrongVerifiedDictionary$wrongOrganismsVerified) |>
  dplyr::arrange(dplyr::desc(sortScore)) |>
  dplyr::distinct(name, dataSourceTitleShort, .keep_all = TRUE) |>
  dplyr::select(
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

dataCleanedOrganismVerified <- dplyr::left_join(
  dataCleanedOrganism,
  verified_df
) |>
  dplyr::select(
    organismType,
    organismValue,
    organismCleaned,
    organismDbTaxo,
    taxonId,
    currentName,
    currentCanonicalFull,
    taxonomy,
    rank
  ) |>
  dplyr::filter(!is.na(organismDbTaxo))

## last version does not contain "species"
indexFungorum <- dataCleanedOrganismVerified |>
  dplyr::filter(organismDbTaxo == "Index Fungorum") |>
  dplyr::mutate(
    rank = ifelse(
      test = rank == "kingdom|phylum|class|order|family|",
      yes = "kingdom|phylum|class|order|family|species",
      no = rank
    )
  )

if (nrow(indexFungorum != 0)) {
  dataCleanedOrganismVerified <- dataCleanedOrganismVerified |>
    dplyr::filter(organismDbTaxo != "Index Fungorum") |>
    dplyr::bind_rows(indexFungorum)
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
  dataCleanedOrganismManipulated <- data.frame() |>
    dplyr::mutate(
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
  dataCleanedOrganismManipulated |>
  dplyr::distinct(organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    .keep_all = TRUE
  ) |>
  dplyr::select(
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    organismCleanedRank
  )

## Avoid generic homonyms
dataCleanedOrganismManipulated_clean <-
  dataCleanedOrganismManipulated_clean |>
  dplyr::mutate(n = stringr::str_count(
    string = organismDetected,
    pattern = stringr::fixed(" ")
  )) |>
  dplyr::filter(n != 0 | organismCleaned == organismDetected) |>
  dplyr::select(-n)

dataCleanedOrganismManipulated_clean_2 <-
  dplyr::left_join(
    dataCleanedOrganismManipulated_clean,
    dataCleanedOrganismManipulated
  ) |>
  dplyr::distinct()

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

dataCuratedOrganismAuto <- dataCuratedOrganismAuto |>
  dplyr::select(
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
  ) |>
  dplyr::filter(grepl(pattern = "[[:alnum:]]", x = organismTaxonRanks))

dataCuratedOrganismAuto <- dataCuratedOrganismAuto |>
  dplyr::anti_join(
    wrongHomonymsDictionary |>
      dplyr::distinct(organismCleaned,
        organismCleanedRank = organismCleaned_rank
      )
  )

log_debug("exporting ... ")
log_debug(pathDataInterimTablesProcessedOrganismFinal)
readr::write_delim(
  x = dataCuratedOrganismAuto,
  delim = "\t",
  file = pathDataInterimTablesProcessedOrganismFinal,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
