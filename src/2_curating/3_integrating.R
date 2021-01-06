cat("This script integrates all results together. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
source("r/vroom_safe.R")

cat("loading files ... \n")
cat("... original table \n")
originalTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalTable)

originalStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesOriginalStructureFull)

cat("loading dictionaries ... \n")
if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
  cat("... structures \n")
  structureDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)
}

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  cat("... organisms \n")
  organismDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  cat("... references \n")
  referenceOrganismDictionary <-
    vroom_read_safe(path = pathDataInterimDictionariesReferenceOrganismDictionary)
}

if (file.exists(pathDataInterimDictionariesStructureMetadata)) {
  cat("... structures metadata \n")
  structureMetadata <-
    vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata)
}

if (file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  cat("... organisms metadata \n")
  organismMetadata <-
    vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata)
}

cat("... cleaned organisms \n")
organismTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedOrganismFinal) %>%
  select(
    organismOriginal,
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

cat("... translated structures \n")
translatedStructureTable <-
  vroom_read_safe(path = pathDataInterimTablesTranslatedStructureFinal)

cat("... cleaned structures \n")
cleanedStructureTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedStructureNamed) %>%
  select(
    structureTranslated,
    structureCleanedSmiles = smilesSanitized,
    structureCleanedInchi = inchiSanitized,
    structureCleanedInchikey3D = inchikeySanitized,
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

cat("... cleaned references \n")
referenceTableFull <-
  vroom_read_safe(path = pathDataInterimTablesCleanedReferenceFile) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

cat("joining ... \n")
if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
  file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  cat("... previously cleaned organisms with metadata \n")
  organismOld <-
    left_join(organismDictionary, organismMetadata)
}

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
  file.exists(pathDataInterimDictionariesStructureMetadata)) {
  cat("... previously cleaned structures with metadata \n")
  structureOld <-
    left_join(structureDictionary, structureMetadata)
}

if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
  file.exists(pathDataInterimDictionariesOrganismMetadata)) {
  cat("... previously cleaned organism with new ones \n")
  organismTableFull <- bind_rows(organismTableFull, organismOld) %>%
    distinct()
}

cat("... translated structures with cleaned ones ... \n")
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated)

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
  file.exists(pathDataInterimDictionariesStructureMetadata)) {
  cat("... previously cleaned structures \n")
  structureFull <- bind_rows(structureFull, structureOld) %>%
    distinct()
}

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary)) {
  cat("... previously cleaned references \n")
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

cat("splitting metadata from minimal columns ... \n")
cat("... structures \n")
structureMinimal <- structureFull %>%
  filter(!is.na(structureCleanedInchikey3D)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles
  )

structureMetadata <- structureFull %>%
  filter(!is.na(structureCleanedInchikey3D)) %>%
  distinct(
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    structureCleaned_inchikey2D,
    structureCleaned_validatorLog,
    structureCleaned_molecularFormula,
    structureCleaned_exactMass,
    structureCleaned_xlogp,
    structureCleaned_stereocenters_unspecified,
    structureCleaned_stereocenters_total,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  )

cat("... organisms \n")
organismMinimal <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  filter(grepl(pattern = "[A-Za-z]", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(
    organismOriginal,
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

cat("... references \n")
referenceMinimal <- referenceTableFull %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  filter(!is.na(referenceCleanedTitle)) %>%
  distinct(
    organismOriginal,
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
    organismOriginal,
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

cat("cleaning memory ... \n")
gc(
  verbose = TRUE,
  reset = TRUE,
  full = TRUE
)
rm(organismTableFull)

cat("joining minimal table ... \n")
inhouseDbMinimal <-
  left_join(originalTable, structureMinimal) %>%
  filter(!is.na(structureCleanedInchikey3D)) %>%
  left_join(., organismMinimal) %>%
  filter(!is.na(organismCleaned)) %>%
  left_join(., referenceMinimal) %>%
  filter(!is.na(referenceCleanedTitle) |
    database == "dnp_1") %>%
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid) |
      database == "dnp_1"
  ) %>%
  distinct(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    # structureCleanedName,
    # structureCleanedNameIupac,
    referenceCleanedTitle
  )

cat("outputting table with missing empty translations (for later on) ... \n")
openDbMaximal <- originalTable %>%
  distinct(
    database,
    organismOriginal,
    referenceType,
    referenceValue,
    structureType,
    structureValue
  ) %>%
  filter(referenceType != "authors") %>%
  filter(referenceType != "journal") %>%
  filter(referenceType != "external") %>%
  filter(referenceType != "isbn")

cat(
  "generating list with chemical names having no translation \n",
  "to avoid translating them again (since process is long) \n"
)
structureNA <- anti_join(
  x = originalTable,
  y = structureMinimal
) %>%
  distinct(
    structureType,
    structureValue
  )

structureNA <- left_join(structureNA, structureFull) %>%
  filter(is.na(structureCleanedInchikey3D)) %>%
  distinct(structureType,
    structureValue,
    .keep_all = TRUE
  ) %>%
  select(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    # structureCleanedName,
    # structureCleanedNameIupac
  ) %>%
  distinct()

cat("ensuring directories exist \n")
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

cat("writing the monster table, if running fullmode, this may take a while \n")
cat(pathDataInterimTablesCuratedTable, "\n")
vroom_write_safe(
  x = inhouseDbMinimal,
  path = pathDataInterimTablesCuratedTable
)

cat(pathDataInterimDictionariesStructureDictionary, "\n")
vroom_write_safe(
  x = structureMinimal,
  path = pathDataInterimDictionariesStructureDictionary
)

cat(pathDataInterimDictionariesStructureAntiDictionary, "\n")
vroom_write_safe(
  x = structureNA,
  path = pathDataInterimDictionariesStructureAntiDictionary
)

cat(pathDataInterimDictionariesOrganismDictionary, "\n")
vroom_write_safe(
  x = organismMinimal,
  path = pathDataInterimDictionariesOrganismDictionary
)

cat(
  pathDataInterimDictionariesReferenceOrganismDictionary,
  "\n"
)
vroom_write_safe(
  x = referenceTableFull,
  path = pathDataInterimDictionariesReferenceOrganismDictionary
)

cat(pathDataInterimDictionariesStructureMetadata, "\n")
vroom_write_safe(
  x = structureMetadata,
  path = pathDataInterimDictionariesStructureMetadata
)

cat(pathDataInterimDictionariesOrganismMetadata, "\n")
vroom_write_safe(
  x = organismMetadata,
  path = pathDataInterimDictionariesOrganismMetadata
)

cat(pathDataInterimDictionariesReferenceMetadata, "\n")
vroom_write_safe(
  x = referenceMetadata,
  path = pathDataInterimDictionariesReferenceMetadata
)

cat(pathDataInterimTablesCuratedTableMaximal, "\n")
vroom_write_safe(
  x = openDbMaximal,
  path = pathDataInterimTablesCuratedTableMaximal
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
