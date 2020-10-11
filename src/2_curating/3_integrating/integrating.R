cat("This script integrates all results together. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/analysis.R")
source("functions/helpers.R")

cat("loading files ... \n")
cat("... original table \n")
originalTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesOriginalTable),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("loading dictionaries ... \n")
if (file.exists(pathDataInterimDictionariesStructureDictionary))
  cat("... structures \n")

if (file.exists(pathDataInterimDictionariesStructureDictionary))
  structureDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesStructureDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

if (file.exists(pathDataInterimDictionariesOrganismDictionary))
  cat("... organisms \n")

if (file.exists(pathDataInterimDictionariesOrganismDictionary))
  organismDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesOrganismDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary))
  cat("... references \n")

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary))
  referenceOrganismDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesReferenceOrganismDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

if (file.exists(pathDataInterimDictionariesStructureMetadata))
  cat("... structures metadata \n")

if (file.exists(pathDataInterimDictionariesStructureMetadata))
  structureMetadata <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesStructureMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

if (file.exists(pathDataInterimDictionariesOrganismMetadata))
  cat("... organisms metadata \n")

if (file.exists(pathDataInterimDictionariesOrganismMetadata))
  organismMetadata <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesOrganismMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

cat("... cleaned organisms \n")
organismTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedOrganismFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    organismOriginal,
    organismCleaned,
    organismCleaned_dbTaxo = organismDbTaxo,
    organismCleaned_dbTaxoTaxonIds = organismTaxonIds,
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
  )

cat("... translated structures \n")
translatedStructureTable <- read_delim(
  file = gzfile(description = pathDataInterimTablesTranslatedStructureFinal),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("... cleaned structures \n")
cleanedStructureTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    structureTranslated,
    structureCleanedSmiles = smilesSanitized,
    structureCleanedInchi = inchiSanitized,
    structureCleanedInchikey3D = inchikeySanitized,
    structureCleaned_inchikey2D = shortikSanitized,
    structureCleaned_validatorLog = validatorLog,
    structureCleaned_molecularFormula = formulaSanitized,
    structureCleaned_exactMass = exactmassSanitized,
    structureCleaned_xlogp = xlogpSanitized
  )

cat("... classified structures \n")
classifiedStructureTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedStructureFileClassified),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    structureCleanedSmiles = smiles,
    structureCleaned_class = class_results,
    structureCleaned_superclass = superclass_results,
    structureCleaned_pathway = pathway_results,
    structureCleaned_glycoside = isglycoside,
    structureCleaned_fp1 = fp1,
    structureCleaned_fp2 = fp2
  )

cat("... cleaned references \n")
referenceTableFull <- read_delim(
  file = gzfile(description = pathDataInterimTablesCleanedReferenceFile),
  delim = "\t",
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("joining ... \n")
if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata))
  cat("... previously cleaned organisms with metadata \n")

if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata))
  organismOld <-
  left_join(organismDictionary, organismMetadata)

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata))
  cat("... previously cleaned structures with metadata \n")

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata))
  structureOld <-
  left_join(structureDictionary, structureMetadata)

if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata))
  cat("... previously cleaned organism with new ones \n")

if (file.exists(pathDataInterimDictionariesOrganismDictionary) &
    file.exists(pathDataInterimDictionariesOrganismMetadata))
  organismTableFull <- bind_rows(organismTableFull, organismOld) %>%
  distinct()

cat("... translated structures with cleaned ones ... \n")
structureFull <-
  left_join(translatedStructureTable, cleanedStructureTableFull) %>%
  select(-structureTranslated)

cat("... with their classification \n")
structureFull <-
  left_join(structureFull, classifiedStructureTableFull)

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata))
  cat("... previously cleaned and classified structures \n")

if (file.exists(pathDataInterimDictionariesStructureDictionary) &
    file.exists(pathDataInterimDictionariesStructureMetadata))
  structureFull <- bind_rows(structureFull, structureOld) %>%
  distinct()

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary))
  cat("... previously cleaned references \n")

if (file.exists(pathDataInterimDictionariesReferenceOrganismDictionary))
  referenceTableFull <-
  bind_rows(referenceTableFull, referenceOrganismDictionary) %>%
  distinct()

cat("splitting metadata from minimal columns ... \n")
cat("... structures \n")
structureMinimal <- structureFull %>%
  filter(!is.na(structureCleanedInchi)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    # structureCleanedName,
    # structureCleanedNameIupac
  )

structureMetadata <- structureFull %>%
  filter(!is.na(structureCleanedInchi)) %>%
  distinct(
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    structureCleaned_inchikey2D,
    structureCleaned_validatorLog,
    structureCleaned_molecularFormula,
    structureCleaned_exactMass,
    structureCleaned_xlogp,
    # structureCleaned_class,
    # structureCleaned_superclass,
    # structureCleaned_pathway,
    # structureCleaned_glycoside,
    # structureCleaned_fp1,
    # structureCleaned_fp2
  )

cat("... organisms \n")
organismMinimal <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  filter(grepl(pattern = "[A-Za-z]", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(
    organismOriginal,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

organismMetadata <- organismTableFull %>%
  filter(!is.na(organismCleaned)) %>%
  filter(grepl(pattern = "[A-Za-z]", x = organismCleaned_dbTaxoTaxonRanks)) %>%
  distinct(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
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
    organismCleaned,
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
    organismCleaned,
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

# cleaning memory
gc(verbose = TRUE,
   reset = TRUE,
   full = TRUE)

cat("joining minimal table ... \n")
cat("... structures \n")
inhouseDbMinimal <-
  left_join(originalTable, structureMinimal) %>%
  filter(!is.na(structureCleanedInchikey3D))

cat("... organisms \n")
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, organismMinimal) %>%
  filter(!is.na(organismCleaned))

cat("... references \n")
inhouseDbMinimal <-
  left_join(inhouseDbMinimal, referenceMinimal) %>%
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
structureNA <- anti_join(x = originalTable,
                         y = structureFull)

structureNA <- left_join(structureNA, structureFull) %>%
  filter(is.na(structureCleanedInchi)) %>%
  distinct(
    structureType,
    structureValue,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    # structureCleanedName,
    # structureCleanedNameIupac
  )

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

cat("writing the monster table, if running fullmode, this may take a while \n")
cat(pathDataInterimTablesCuratedTable, "\n")
write.table(
  x = inhouseDbMinimal,
  file = gzfile(
    description = pathDataInterimTablesCuratedTable,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesStructureDictionary, "\n")
write.table(
  x = structureMinimal,
  file = gzfile(
    description = pathDataInterimDictionariesStructureDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesStructureAntiDictionary, "\n")
write.table(
  x = structureNA,
  file = gzfile(
    description = pathDataInterimDictionariesStructureAntiDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesOrganismDictionary, "\n")
write.table(
  x = organismMinimal,
  file = gzfile(
    description = pathDataInterimDictionariesOrganismDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesReferenceOrganismDictionary,
    "\n")
write.table(
  x = referenceTableFull,
  file = gzfile(
    description = pathDataInterimDictionariesReferenceOrganismDictionary,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesStructureMetadata, "\n")
write.table(
  x = structureMetadata,
  file = gzfile(
    description = pathDataInterimDictionariesStructureMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesOrganismMetadata, "\n")
write.table(
  x = organismMetadata,
  file = gzfile(
    description = pathDataInterimDictionariesOrganismMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimDictionariesReferenceMetadata, "\n")
write.table(
  x = referenceMetadata,
  file = gzfile(
    description = pathDataInterimDictionariesReferenceMetadata,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesCuratedTableMaximal, "\n")
write.table(
  x = openDbMaximal,
  file = gzfile(
    description = pathDataInterimTablesCuratedTableMaximal,
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
