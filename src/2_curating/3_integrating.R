source("r/log_debug.R")
log_debug("This script integrates all results together.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(tidyverse)
source("r/vroom_safe.R")

log_debug("loading files ...")
log_debug("... original table")
originalTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalTable)

originalStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalStructureFull)

log_debug("loading dictionaries ...")
if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
  log_debug("... structures")
  structureDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)
}

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  log_debug("... organisms")
  organismDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  log_debug("... references")
  referenceOrganismDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesReferenceOrganismDictionary)
}

if (file.exists(pathDataInterimDictionariesStructureMetadata)) {
  log_debug("... structures metadata")
  structureMetadata <-
    vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata)
}

if (file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  log_debug("... organisms metadata")
  organismMetadata <-
    vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata)
}

log_debug("... cleaned organisms")
organismTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismFinal) %>%
  select(
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
  ) %>%
  distinct()

log_debug("... translated structures")
translatedStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedStructureFinal)

log_debug("... cleaned structures")
cleanedStructureTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureNamed) %>%
  select(
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
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  )

log_debug("... cleaned references")
referenceTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedReferenceFile) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

log_debug("joining ...")
if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
  file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  log_debug("... previously cleaned organisms with metadata")
  organismOld <-
    left_join(organismDictionary, organismMetadata)
}

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
  file.exists(pathDataInterimDictionariesStructureMetadata)) {
  log_debug("... previously cleaned structures with metadata")
  structureOld <-
    left_join(structureDictionary, structureMetadata)
}

if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
  file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  log_debug("... previously cleaned organism with new ones")
  organismTableFull <- bind_rows(organismTableFull, organismOld) %>%
    distinct()
}

log_debug("... translated structures with cleaned ones ...")
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated) %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey,
    .keep_all = TRUE
  )

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
  file.exists(pathDataInterimDictionariesStructureMetadata)) {
  log_debug("... previously cleaned structures")
  structureFull <- bind_rows(structureFull, structureOld) %>%
    distinct(
      structureType,
      structureValue,
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey,
      .keep_all = TRUE
    )
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  log_debug("... previously cleaned references")
  referenceTableFull <-
    bind_rows(referenceTableFull, referenceOrganismDictionary) %>%
    distinct()
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
structureMinimal <- structureFull %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles
  )

structureMetadata <- structureFull %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  select(-structureType, -structureValue) %>%
  distinct(
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    .keep_all = TRUE
  )

log_debug("... organisms")
organismMinimal <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  filter(grepl(pattern = "[A-Za-z]", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(
    organismType,
    organismValue,
    organismDetected,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_id,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

organismMetadata <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  filter(grepl(pattern = "[A-Za-z]", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(
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
referenceMinimal <- referenceTableFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  filter(!is.na(referenceCleanedTitle)) %>%
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

referenceMetadata <- referenceTableFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  filter(!is.na(referenceCleanedTitle)) %>%
  distinct(
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

log_debug("cleaning memory ...")
gc(
  verbose = TRUE,
  reset = TRUE,
  full = TRUE
)
rm(organismTableFull)

log_debug("joining minimal table ...")
inhouseDbMinimal <-
  left_join(originalTable, structureMinimal) %>%
  filter(!is.na(structureCleanedInchikey)) %>%
  left_join(., organismMinimal) %>%
  filter(!is.na(organismCleaned)) %>%
  left_join(., referenceMinimal %>%
    select(-organismDetected)) %>%
  filter(!is.na(referenceCleanedTitle) |
    database == "dnp") %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid) |
      database == "dnp"
  ) %>%
  distinct(
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
    # structureCleanedName,
    # structureCleanedNameIupac,
    referenceCleanedTitle
  )

log_debug("outputting table with missing empty translations (for later on) ...")
openDbMaximal <- originalTable %>%
  distinct(
    database,
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    structureType,
    structureValue
  ) %>%
  filter(referenceType != "authors") %>%
  filter(referenceType != "journal") %>%
  filter(referenceType != "external") %>%
  filter(referenceType != "isbn")

log_debug(
  "generating list with chemical names having no translation \n",
  "to avoid translating them again (since process is long)"
)
structureNA <- anti_join(
  x = originalStructureTable,
  y = structureMinimal
) %>%
  distinct(
    structureType,
    structureValue
  )

structureNA <- left_join(structureNA, structureFull) %>%
  filter(is.na(structureCleanedInchikey)) %>%
  distinct(structureType,
    structureValue,
    .keep_all = TRUE
  ) %>%
  select(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    # structureCleanedName,
    # structureCleanedNameIupac
  ) %>%
  distinct()

log_debug("ensuring directories exist")
ifelse(
  test = !dir.exists(pathDataInterimTablesCurated),
  yes = dir.create(pathDataInterimTablesCurated),
  no = paste(pathDataInterimTablesCurated, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimDictionariesStructure),
  yes = dir.create(pathDataInterimDictionariesStructure),
  no = paste(pathDataInterimDictionariesStructure, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimDictionariesOrganism),
  yes = dir.create(pathDataInterimDictionariesOrganism),
  no = paste(pathDataInterimDictionariesOrganism, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimDictionariesReference),
  yes = dir.create(pathDataInterimDictionariesReference),
  no = paste(pathDataInterimDictionariesReference, "exists")
)

log_debug("writing the monster table, if running fullmode, this may take a while")
log_debug(pathDataInterimTablesCuratedTable)
vroom_write_safe_append(
  x = inhouseDbMinimal,
  path = pathDataInterimTablesCuratedTable
)

log_debug(pathDataInterimDictionariesStructureDictionary)
vroom_write_safe(
  x = structureMinimal,
  path = pathDataInterimDictionariesStructureDictionary
)

log_debug(pathDataInterimDictionariesStructureAntiDictionary)
vroom_write_safe_append(
  x = structureNA,
  path = pathDataInterimDictionariesStructureAntiDictionary
)

log_debug(pathDataInterimDictionariesOrganismDictionary)
vroom_write_safe(
  x = organismMinimal,
  path = pathDataInterimDictionariesOrganismDictionary
)

log_debug(
  pathDataInterimDictionariesReferenceOrganismDictionary
)
vroom_write_safe(
  x = referenceTableFull,
  path = pathDataInterimDictionariesReferenceOrganismDictionary
)

log_debug(pathDataInterimDictionariesStructureMetadata)
vroom_write_safe(
  x = structureMetadata,
  path = pathDataInterimDictionariesStructureMetadata
)

log_debug(pathDataInterimDictionariesOrganismMetadata)
vroom_write_safe(
  x = organismMetadata,
  path = pathDataInterimDictionariesOrganismMetadata
)

log_debug(pathDataInterimDictionariesReferenceMetadata)
vroom_write_safe(
  x = referenceMetadata,
  path = pathDataInterimDictionariesReferenceMetadata
)

log_debug(pathDataInterimTablesCuratedTableMaximal)
vroom_write_safe_append(
  x = openDbMaximal,
  path = pathDataInterimTablesCuratedTableMaximal
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
