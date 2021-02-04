cat("This script is the first attempt to create the tables for sql use \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(data.table)
library(DBI)
library(RSQLite)
library(tidyverse)
source("r/vroom_safe.R")
source("r/sqlFromFile.R")
source("r/dbSendQueries.R")

## chemical support tables to discuss and do
## references support tables to discuss and do

dbTypes <- read_delim(
  file = "../docs/dataset.tsv",
  delim = "\t"
) %>%
  select(
    database,
    type
  )

dbList <- lapply(pathDataInterimDbDir, vroom_read_safe)

structureDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)

organismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

referenceOrganismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesReferenceOrganismDictionary)

structureMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata)

organismMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata)

inhouseDbMinimal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTable)

manuallyValidated <-
  vroom_read_safe(path = "../data/validation/manuallyValidated.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_1 = 1)

automaticallyValidated <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_2 = 2)

manuallyRemoved <-
  vroom_read_safe(path = "../data/validation/manuallyRemoved.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationStateId_3 = 3)

curation_states <-
  data.frame(
    id = c(1, 2, 3, 4),
    name = c(
      "manual_validation",
      "automatic_validation",
      "manual_removal",
      "automatic_removal"
    )
  )

structureOld <-
  left_join(
    structureDictionary,
    structureMetadata
  )

organismOld <-
  left_join(
    organismDictionary,
    organismMetadata
  )

dbTable <- rbindlist(l = dbList, fill = TRUE) %>%
  select(
    database,
    organismOriginal_clean = organism_clean,
    organismOriginal_dirty = organism_dirty,
    structureOriginal_inchi = structure_inchi,
    structureOriginal_nominal = structure_name,
    structureOriginal_smiles = structure_smiles,
    referenceOriginal_authors = reference_authors,
    referenceOriginal_doi = reference_doi,
    referenceOriginal_external = reference_external,
    referenceOriginal_isbn = reference_isbn,
    referenceOriginal_journal = reference_journal,
    referenceOriginal_original = reference_original,
    referenceOriginal_pubmed = reference_pubmed,
    referenceOriginal_publishingDetails = reference_publishingDetails,
    referenceOriginal_split = reference_split,
    referenceOriginal_title = reference_title,
  ) %>%
  tibble()

originalTable <- dbTable %>%
  pivot_longer(
    cols = 7:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 4:6,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 2:3,
    names_to = c("drop3", "organismType"),
    names_sep = "_",
    values_to = "organismValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop, -drop2, -drop3) %>%
  distinct()

database_source <- originalTable %>%
  distinct(database) %>%
  mutate(databaseSourceId = row_number()) %>%
  left_join(., dbTypes) %>%
  group_by(type) %>%
  mutate(typeId = cur_group_id()) %>%
  ungroup()

organism_source <- originalTable %>%
  distinct(organismType, organismValue) %>%
  mutate(organismSourceId = row_number()) %>%
  group_by(organismType) %>%
  mutate(typeId = cur_group_id()) %>%
  ungroup()

structure_source <- originalTable %>%
  distinct(structureType, structureValue) %>%
  mutate(structureSourceId = row_number()) %>%
  group_by(structureType) %>%
  mutate(typeId = cur_group_id()) %>%
  ungroup()

reference_source <- originalTable %>%
  distinct(referenceType, referenceValue) %>%
  mutate(referenceSourceId = row_number()) %>%
  group_by(referenceType) %>%
  mutate(typeId = cur_group_id()) %>%
  ungroup()

database_types <- database_source %>%
  select(
    id = typeId,
    name = type
  ) %>%
  distinct() %>%
  arrange(id)

organism_types <- organism_source %>%
  select(
    id = typeId,
    name = organismType
  ) %>%
  distinct() %>%
  arrange(id)

structure_types <- structure_source %>%
  select(
    id = typeId,
    name = structureType
  ) %>%
  distinct() %>%
  arrange(id)

reference_types <- reference_source %>%
  select(
    id = typeId,
    name = referenceType
  ) %>%
  distinct() %>%
  arrange(id)

data_source <- originalTable %>%
  left_join(., database_source %>% select(-typeId)) %>%
  left_join(., organism_source %>% select(-typeId)) %>%
  left_join(., structure_source %>% select(-typeId)) %>%
  left_join(., reference_source %>% select(-typeId)) %>%
  mutate(id = row_number())

database_source <- database_source %>%
  select(
    id = databaseSourceId,
    name = database,
    typeId
  )

organism_source <- organism_source %>%
  select(
    id = organismSourceId,
    value = organismValue,
    typeId
  )

structure_source <- structure_source %>%
  select(
    id = structureSourceId,
    value = structureValue,
    typeId
  )

reference_source <- reference_source %>%
  select(
    id = referenceSourceId,
    value = referenceValue,
    typeId
  )

organism_detected <- organismOld %>%
  distinct(
    organismDetected,
    organismCleaned
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    organismDetected,
    organismCleaned
  )

organism_cleaned <- organismOld %>%
  distinct(organismCleaned) %>%
  mutate(id = row_number()) %>%
  select(id,
    name = organismCleaned
  )

organism_synonyms <- organism_detected %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    id = id.x,
    name = organismDetected,
    organismCleanedId = id.y
  )

taxonomic_databases <- organismOld %>%
  distinct(organismCleaned_dbTaxo) %>%
  group_by(organismCleaned_dbTaxo) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  select(id,
    name = organismCleaned_dbTaxo
  )

taxonomic_information <- organismOld %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    cleanedOrganismId = id,
    everything()
  ) %>%
  left_join(.,
    taxonomic_databases,
    by = c("organismCleaned_dbTaxo" = "name")
  ) %>%
  select(
    taxonomicDatabaseId = id,
    everything()
  ) %>%
  distinct(cleanedOrganismId,
    taxonomicDatabaseId,
    .keep_all = TRUE
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    cleanedOrganismId,
    taxonomicDatabaseId,
    taxonomicId = organismCleaned_id,
    ranks = organismCleaned_dbTaxoTaxonRanks,
    taxonomy = organismCleaned_dbTaxoTaxonomy,
    rank = organismCleaned_rank
  )

reference_cleaned <- referenceOrganismDictionary %>%
  distinct(
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  ) %>%
  mutate(id = row_number()) %>%
  left_join(., referenceOrganismDictionary) %>%
  select(
    id,
    doi = referenceCleanedDoi,
    pmcid = referenceCleanedPmcid,
    pmid = referenceCleanedPmid,
    title = referenceCleanedTitle,
    everything()
  )

structure_cleaned <- structureOld %>%
  distinct(
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey
  ) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey
  ) %>%
  left_join(
    .,
    structureOld %>% distinct(
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey,
      .keep_all = TRUE
    )
  ) %>%
  distinct() %>%
  select(
    id,
    traditionalName = structureCleaned_nameTraditional,
    iupacName = structureCleaned_nameIupac,
    inchikey = structureCleanedInchikey,
    inchikey2D = structureCleaned_inchikey2D,
    inchi = structureCleanedInchi,
    inchi2D = structureCleaned_inchi2D,
    smiles = structureCleanedSmiles,
    smiles2D = structureCleaned_smiles2D,
    stereocentersTotal = structureCleaned_stereocenters_total,
    stereocentersUnspecified = structureCleaned_stereocenters_unspecified,
    molecularFormula = structureCleaned_molecularFormula,
    exactMass = structureCleaned_exactMass,
    xlogp = structureCleaned_xlogp
  )

inhouseDbMinimal_complemented <- inhouseDbMinimal %>%
  mutate(curationStateId_4 = 4) %>%
  left_join(., manuallyValidated) %>%
  left_join(., manuallyRemoved) %>%
  left_join(., automaticallyValidated) %>%
  pivot_longer(cols = (ncol(.) - 3):ncol(.)) %>%
  arrange(value) %>%
  distinct(
    database,
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    structureType,
    structureValue,
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey,
    organismCleaned,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  select(everything(),
    -name,
    curationStateId = value
  )

data_processed_temp <- inhouseDbMinimal_complemented %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    organismCleanedId = id,
    everything()
  ) %>%
  left_join(
    .,
    structure_cleaned %>%
      distinct(
        id,
        inchikey,
        inchi,
        smiles
      ),
    by = c(
      "structureCleanedSmiles" = "smiles",
      "structureCleanedInchi" = "inchi",
      "structureCleanedInchikey" = "inchikey"
    )
  ) %>%
  select(
    structureCleanedId = id,
    everything()
  ) %>%
  left_join(
    .,
    reference_cleaned,
    by = c(
      "organismType" = "organismType",
      "organismValue" = "organismValue",
      "referenceType" = "referenceType",
      "referenceValue" = "referenceValue",
      "referenceCleanedTitle" = "title"
    )
  ) %>%
  select(
    referenceCleanedId = id,
    everything()
  ) %>%
  distinct(
    database,
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    structureType,
    structureValue,
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId,
    curationStateId
  )

data_processed <- data_processed_temp %>%
  distinct(
    organismCleanedId,
    structureCleanedId,
    referenceCleanedId,
    curationStateId
  ) %>%
  mutate(id = row_number()) %>%
  left_join(., data_processed_temp)

data_processed__data_source <- data_processed_temp %>%
  left_join(., data_source %>%
    select(
      dataSourceId = id,
      everything(),
      -type
    )) %>%
  left_join(., data_processed %>%
    select(
      dataProcessedId = id,
      everything()
    )) %>%
  filter(!is.na(dataProcessedId)) %>%
  distinct(dataSourceId, dataProcessedId, curationStateId) %>%
  mutate(id = row_number()) %>%
  select(
    id,
    dataSourceId,
    dataProcessedId,
    curationStateId
  )

data_processed <- data_processed %>%
  distinct(
    id,
    structureCleanedId,
    organismCleanedId,
    referenceCleanedId
  ) %>%
  select(
    id,
    structureCleanedId,
    organismCleanedId,
    referenceCleanedId
  )

data_source <- data_source %>%
  select(
    id,
    databaseSourceId,
    organismSourceId,
    structureSourceId,
    referenceSourceId
  )

reference_cleaned <- reference_cleaned %>%
  distinct(
    id,
    doi,
    pmcid,
    pmid,
    title
  )

rm(
  automaticallyValidated,
  dbList,
  dbTable,
  dbTypes,
  inhouseDbMinimal,
  inhouseDbMinimal_complemented,
  manuallyRemoved,
  manuallyValidated,
  organismDictionary,
  organismMetadata,
  organism_detected,
  organismOld,
  originalTable,
  referenceOrganismDictionary,
  structureDictionary,
  structureMetadata,
  structureOld
)

drv <- SQLite()

## TEMP
file.create(lotusDB)

## TEMP
db <- dbConnect(
  drv = drv,
  dbname = lotusDB
)
## TEMP
dbSendQueries(
  conn = db,
  sql = sqlFromFile("schema_db/0000_create_initial_tables.sql")
)

dbListObjects(db)

dbListFields(db, "chemical_databases")
# colnames(chemical_databases) ## TO DO
dbListFields(db, "chemical_information")
# colnames(chemical_information) ## TO DO
dbListFields(db, "curation_states")
colnames(curation_states)
dbListFields(db, "data_processed")
colnames(data_processed)
dbListFields(db, "data_processed__data_source")
colnames(data_processed__data_source)
dbListFields(db, "data_source")
colnames(data_source)
dbListFields(db, "database_source")
colnames(database_source)
dbListFields(db, "database_types")
colnames(database_types)
dbListFields(db, "organism_cleaned")
colnames(organism_cleaned)
dbListFields(db, "organism_source")
colnames(organism_source)
dbListFields(db, "organism_synonyms")
colnames(organism_synonyms)
dbListFields(db, "organism_types")
colnames(organism_types)
dbListFields(db, "reference_cleaned")
colnames(reference_cleaned)
dbListFields(db, "reference_databases")
# colnames(reference_databases) ## TO DO
dbListFields(db, "reference_information")
# colnames(reference_information) ## TO DO
colnames(reference_cleaned)
dbListFields(db, "reference_source")
colnames(reference_source)
dbListFields(db, "reference_types")
colnames(reference_types)
dbListFields(db, "structure_cleaned")
colnames(structure_cleaned)
dbListFields(db, "structure_source")
colnames(structure_source)
dbListFields(db, "structure_types")
colnames(structure_types)
dbListFields(db, "taxonomic_databases")
colnames(taxonomic_databases)
dbListFields(db, "taxonomic_information")
colnames(taxonomic_information)

# dbWriteTable(
#   conn = db,
#   name = "chemical_databases",
#   value = chemical_databases,
#   row.names = FALSE,
#   append = TRUE
# )

# dbWriteTable(
#   conn = db,
#   name = "chemical_information",
#   value = chemical_information,
#   row.names = FALSE,
#   append = TRUE
# )

dbWriteTable(
  conn = db,
  name = "curation_states",
  value = curation_states,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_processed",
  value = data_processed,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_processed__data_source",
  value = data_processed__data_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "data_source",
  value = data_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "database_source",
  value = database_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "database_types",
  value = database_types,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_cleaned",
  value = organism_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_source",
  value = organism_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_synonyms",
  value = organism_synonyms,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "organism_types",
  value = organism_types,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "reference_cleaned",
  value = reference_cleaned,
  row.names = FALSE,
  append = TRUE
)

# dbWriteTable(
#   conn = db,
#   name = "reference_databases",
#   value = reference_databases,
#   row.names = FALSE,
#   append = TRUE
# )

# dbWriteTable(
#   conn = db,
#   name = "reference_information",
#   value = reference_information,
#   row.names = FALSE,
#   append = TRUE
# )

dbWriteTable(
  conn = db,
  name = "reference_source",
  value = reference_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "reference_types",
  value = reference_types,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structure_cleaned",
  value = structure_cleaned,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structure_source",
  value = structure_source,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "structure_types",
  value = structure_types,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "taxonomic_databases",
  value = taxonomic_databases,
  row.names = FALSE,
  append = TRUE
)

dbWriteTable(
  conn = db,
  name = "taxonomic_information",
  value = taxonomic_information,
  row.names = FALSE,
  append = TRUE
)

dbDisconnect(db)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
