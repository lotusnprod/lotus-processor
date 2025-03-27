source("r/log_debug.R")
log_debug("This script integrates all results together.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)
library(stringr)

log_debug("loading files ...")
log_debug("... original table")
originalTable <-
  readr::read_delim(
    file = pathDataInterimTablesOriginalTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

originalStructureTable <- originalTable |>
  dplyr::distinct(
    structureType,
    structureValue
  )

log_debug("loading dictionaries ...")
if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
  log_debug("... structures")
  structureDictionary <-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureDictionary,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )
}

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("... organisms")
  organismDictionary <-
    readr::read_delim(
      file = pathDataInterimDictionariesOrganismDictionary,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  log_debug("... references")
  referenceOrganismDictionary <-
    readr::read_delim(
      file = pathDataInterimDictionariesReferenceOrganismDictionary,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )
}

if (file.exists(pathDataInterimDictionariesStructureMetadata)) {
  log_debug("... structures metadata")
  structureMetadata <-
    readr::read_delim(
      file = pathDataInterimDictionariesStructureMetadata,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )
}

if (file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  log_debug("... organisms metadata")
  organismMetadata <-
    readr::read_delim(
      file = pathDataInterimDictionariesOrganismMetadata,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )
}

log_debug("... cleaned organisms")
organismTableFull <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedOrganismFinal,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::select(
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    organismCleaned_id = organismCleanedId,
    organismCleaned_rank = organismCleanedRank,
    organismCleaned_dbTaxo = organismDbTaxo,
    # organismCleaned_dbTaxoTaxonIds = organismTaxonIds,
    organismCleaned_dbTaxoTaxonRanks = organismTaxonRanks,
    organismCleaned_dbTaxoTaxonomy = organismTaxonomy,
    organismCleaned_dbTaxo_1kingdom = organism_1_kingdom,
    organismCleaned_dbTaxo_2phylum = organism_2_phylum,
    organismCleaned_dbTaxo_3class = organism_3_class,
    organismCleaned_dbTaxo_4order = organism_4_order,
    organismCleaned_dbTaxo_5family = organism_5_family,
    organismCleaned_dbTaxo_6genus = organism_6_genus,
    organismCleaned_dbTaxo_7species = organism_7_species,
    organismCleaned_dbTaxo_8variety = organism_8_variety
  ) |>
  dplyr::distinct()

log_debug("... translated structures")
translatedStructureTable <-
  readr::read_delim(
    file = pathDataInterimTablesTranslatedStructureFinal,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... cleaned structures")
cleanedStructureTableFull <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedStructureNamed,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::select(
    structureTranslated,
    structureCleanedSmiles = smilesSanitized,
    structureCleanedInchi = inchiSanitized,
    structureCleanedInchikey = inchikeySanitized,
    structureCleaned_smiles2D = smilesSanitizedFlat,
    structureCleaned_inchi2D = inchiSanitizedFlat,
    structureCleaned_inchikey2D = shortikSanitized,
    structureCleaned_validatorLog = validatorLog,
    structureCleaned_molecularFormula = formulaSanitized,
    structureCleaned_exactMass = exactmassSanitized,
    structureCleaned_xlogp = xlogpSanitized,
    structureCleaned_stereocenters_unspecified = count_unspecified_atomic_stereocenters,
    structureCleaned_stereocenters_total = count_atomic_stereocenters,
    structureCleaned_cid,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  )

log_debug("... cleaned references")
referenceTableFull <-
  readr::read_delim(
    file = pathDataInterimTablesProcessedReferenceFile,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimTablesCurated)
create_dir(export = pathDataInterimDictionariesOrganism)
create_dir(export = pathDataInterimDictionariesReference)
create_dir(export = pathDataInterimDictionariesStructure)

log_debug("joining ...")
if (
  file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata)
) {
  log_debug("... previously cleaned organisms with metadata")
  organismOld <-
    dplyr::left_join(organismDictionary, organismMetadata)
}

if (
  file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata)
) {
  log_debug("... previously cleaned structures with metadata")
  structureOld <-
    dplyr::left_join(structureDictionary, structureMetadata)
}

if (
  file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata)
) {
  log_debug("... previously cleaned organism with new ones")
  organismTableFull <-
    dplyr::bind_rows(organismTableFull, organismOld) |>
    dplyr::distinct()
}

log_debug("... translated structures with cleaned ones ...")
structureFull <-
  dplyr::left_join(translatedStructureTable, cleanedStructureTableFull) |>
  dplyr::select(-structureTranslated) |>
  dplyr::filter(!is.na(structureCleanedInchikey)) |>
  dplyr::distinct(
    structureType,
    structureValue,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleaned_cid,
    .keep_all = TRUE
  )

if (
  file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata)
) {
  log_debug("... previously cleaned structures")
  structureFull <-
    dplyr::bind_rows(structureFull, structureOld) |>
    dplyr::distinct(
      structureType,
      structureValue,
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey,
      structureCleaned_cid,
      .keep_all = TRUE
    )
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  log_debug("... previously cleaned references")
  referenceTableFull <-
    dplyr::bind_rows(referenceTableFull, referenceOrganismDictionary) |>
    dplyr::distinct()
}

rm(
  structureOld,
  structureDictionary,
  organismOld,
  organismDictionary,
  referenceOrganismDictionary
)

log_debug("splitting metadata from minimal columns ...")
log_debug("... structures")
structureMinimal <- structureFull |>
  dplyr::filter(!is.na(structureCleanedInchikey)) |>
  dplyr::distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles
  )

structureMetadata <- structureFull |>
  dplyr::filter(!is.na(structureCleanedInchikey)) |>
  dplyr::select(-structureType, -structureValue) |>
  dplyr::distinct(
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    structureCleaned_cid,
    .keep_all = TRUE
  )

log_debug("... organisms")
organismMinimal <- organismTableFull |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::filter(grepl(
    pattern = "[A-Za-z]",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::distinct(
    organismType,
    organismValue,
    organismCleaned,
    organismDetected
  )

organismMetadata <- organismTableFull |>
  dplyr::filter(!is.na(organismCleaned)) |>
  dplyr::filter(grepl(
    pattern = "[A-Za-z]",
    x = organismCleaned_dbTaxoTaxonRanks
  )) |>
  dplyr::distinct(
    organismCleaned,
    organismCleaned_id,
    organismCleaned_rank,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    organismCleaned_dbTaxo_1kingdom,
    organismCleaned_dbTaxo_2phylum,
    organismCleaned_dbTaxo_3class,
    organismCleaned_dbTaxo_4order,
    organismCleaned_dbTaxo_5family,
    organismCleaned_dbTaxo_6genus,
    organismCleaned_dbTaxo_7species,
    organismCleaned_dbTaxo_8variety
  )

log_debug("... references")
referenceMinimal <- referenceTableFull |>
  dplyr::filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) |>
  filter(!is.na(referenceCleanedTitle)) |>
  distinct(
    organismType,
    organismValue,
    organismDetected,
    referenceType,
    referenceValue,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

referenceMetadata <- referenceTableFull |>
  dplyr::filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) |>
  dplyr::filter(!is.na(referenceCleanedTitle)) |>
  dplyr::distinct(
    organismType,
    organismValue,
    organismDetected,
    referenceType,
    referenceValue,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref,
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism,
    referenceCleaned_score_complementDate,
    referenceCleaned_score_complementAuthor,
    referenceCleaned_score_complementJournal,
    referenceCleaned_score_complementTotal
  )

log_debug(
  "generating list with chemical names having no translation \n",
  "to avoid translating them again (since process is long)"
)
structureNA <- dplyr::anti_join(
  x = originalStructureTable,
  y = structureMinimal
) |>
  dplyr::distinct(
    structureType,
    structureValue
  )

structureNA <- dplyr::left_join(structureNA, structureFull) |>
  dplyr::filter(is.na(structureCleanedInchikey)) |>
  dplyr::distinct(structureType, structureValue, .keep_all = TRUE) |>
  dplyr::select(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    # structureCleanedCid,
    # structureCleanedName,
    # structureCleanedNameIupac
  ) |>
  dplyr::distinct()

log_debug("partial export ...")
log_debug(pathDataInterimDictionariesReferenceOrganismDictionary)
readr::write_delim(
  x = referenceTableFull,
  file = pathDataInterimDictionariesReferenceOrganismDictionary,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesStructureDictionary)
readr::write_delim(
  x = structureMinimal,
  file = pathDataInterimDictionariesStructureDictionary,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesStructureAntiDictionary)
readr::write_delim(
  x = structureNA,
  file = pathDataInterimDictionariesStructureAntiDictionary,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesOrganismDictionary)
readr::write_delim(
  x = organismMinimal,
  file = pathDataInterimDictionariesOrganismDictionary,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesStructureMetadata)
readr::write_delim(
  x = structureMetadata,
  file = pathDataInterimDictionariesStructureMetadata,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesOrganismMetadata)
readr::write_delim(
  x = organismMetadata,
  file = pathDataInterimDictionariesOrganismMetadata,
  delim = "\t",
  na = ""
)

log_debug(pathDataInterimDictionariesReferenceMetadata)
readr::write_delim(
  x = referenceMetadata,
  file = pathDataInterimDictionariesReferenceMetadata,
  delim = "\t",
  na = ""
)

log_debug("cleaning memory ...")
rm(
  organismTableFull,
  structureFull,
  referenceTableFull,
  organismMetadata,
  referenceMetadata,
  structureNA
)
gc(
  verbose = TRUE,
  reset = TRUE,
  full = TRUE
)

log_debug("outputting table with missing empty translations (for later on) ...")
openDbMaximal <- originalTable |>
  dplyr::distinct(
    database,
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    structureType,
    structureValue
  ) |>
  dplyr::filter(referenceType != "authors") |>
  dplyr::filter(referenceType != "journal") |>
  dplyr::filter(referenceType != "external") |>
  dplyr::filter(referenceType != "isbn")

log_debug("partial export ...")
log_debug(pathDataInterimTablesCuratedTableMaximal)
readr::write_delim(
  x = openDbMaximal,
  file = pathDataInterimTablesCuratedTableMaximal,
  delim = "\t",
  na = ""
)

log_debug("cleaning memory ...")
rm(openDbMaximal)
gc(
  verbose = TRUE,
  reset = TRUE,
  full = TRUE
)

log_debug("joining minimal table ...")
inhouseDbMinimal <-
  dplyr::inner_join(originalTable, structureMinimal) |>
  dplyr::inner_join(organismMinimal) |>
  dplyr::left_join(
    referenceMinimal |>
      dplyr::select(-organismDetected)
  )

rm(originalTable)

inhouseDbMinimal_1 <- inhouseDbMinimal |>
  dplyr::filter(database %in% forbidden_export) |>
  dplyr::select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedTitle,
    referenceCleanedDoi,
  ) |>
  dplyr::distinct()

inhouseDbMinimal_2 <- inhouseDbMinimal |>
  dplyr::filter(!is.na(referenceCleanedTitle))

rm(inhouseDbMinimal)

inhouseDbMinimal_2_1 <- inhouseDbMinimal_2 |>
  dplyr::filter(!is.na(referenceCleanedDoi)) |>
  dplyr::select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedTitle,
    referenceCleanedDoi,
  ) |>
  dplyr::distinct()
inhouseDbMinimal_2_2 <- inhouseDbMinimal_2 |>
  dplyr::filter(!is.na(referenceCleanedPmcid)) |>
  dplyr::select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedTitle,
    referenceCleanedDoi,
  ) |>
  dplyr::distinct()
inhouseDbMinimal_2_3 <- inhouseDbMinimal_2 |>
  dplyr::filter(!is.na(referenceCleanedPmid)) |>
  dplyr::select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedTitle,
    referenceCleanedDoi,
  ) |>
  dplyr::distinct()

inhouseDbMinimal <-
  dplyr::bind_rows(
    inhouseDbMinimal_2_1,
    inhouseDbMinimal_2_2,
    inhouseDbMinimal_2_3,
    inhouseDbMinimal_1
  ) |>
  dplyr::distinct()

rm(
  inhouseDbMinimal_2_1,
  inhouseDbMinimal_2_2,
  inhouseDbMinimal_2_3,
  inhouseDbMinimal_1,
  inhouseDbMinimal_2
)

log_debug("removing redundant structures (lower defined stereo)")
inhouseDbMinimal <- inhouseDbMinimal |>
  dplyr::left_join(
    structureMetadata |>
      dplyr::distinct(
        structureCleanedSmiles,
        structureCleanedInchi,
        structureCleanedInchikey,
        structureCleaned_inchikey2D,
        structureCleaned_stereocenters_unspecified,
        structureCleaned_stereocenters_total
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        specified_stereo = as.numeric(structureCleaned_stereocenters_total) -
          as.numeric(structureCleaned_stereocenters_unspecified)
      ) |>
      dplyr::ungroup()
  ) |>
  dplyr::group_by(
    organismCleaned,
    structureCleaned_inchikey2D,
    referenceCleanedTitle,
    referenceCleanedDoi
  ) |>
  dplyr::mutate(best_stereo = max(specified_stereo)) |>
  dplyr::ungroup() |>
  dplyr::filter(specified_stereo == best_stereo) |>
  dplyr::select(
    -structureCleaned_inchikey2D,
    -structureCleaned_stereocenters_unspecified,
    -structureCleaned_stereocenters_total,
    -specified_stereo,
    -best_stereo
  )

log_debug("removing redundant upper taxa")
inhouseDbMinimal <- inhouseDbMinimal |>
  dplyr::mutate(
    temp_org = gsub(
      pattern = " .*",
      replacement = "",
      x = organismCleaned
    ),
    spaces = stringr::str_count(string = organismCleaned, pattern = " ")
  ) |>
  dplyr::mutate(
    spaces = dplyr::if_else(condition = spaces == 0, true = spaces, false = 1)
  ) |>
  ## remove genera where species are found but not species with var
  dplyr::group_by(
    temp_org,
    structureCleanedInchikey,
    referenceCleanedTitle,
    referenceCleanedDoi
  ) |>
  dplyr::mutate(best_org = max(spaces)) |>
  dplyr::ungroup() |>
  dplyr::filter(spaces == best_org) |>
  dplyr::select(-temp_org, -spaces, -best_org)

log_debug(
  "writing the monster table, if running fullmode, this may take a while"
)
log_debug(pathDataInterimTablesCuratedTable)
readr::write_delim(
  x = inhouseDbMinimal,
  file = pathDataInterimTablesCuratedTable,
  delim = "\t",
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
