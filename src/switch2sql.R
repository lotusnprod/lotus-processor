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

if (db_mode == "fromScratch") {
  cat("... creating  the database \n")
  file.create(lotusDB)
}

drv <- SQLite()

cat("... connecting to the database \n")
db <- dbConnect(
  drv = drv,
  dbname = lotusDB
)

if (db_mode == "fromScratch") {
  dbSendQueries(
    conn = db,
    sql = sqlFromFile("schema_db/0000_create_initial_tables.sql")
  )
}

cat("... listing objects \n")
dbListObjects(db)

cat("... loading files ... \n")

dbTypes <- read_delim(
  file = "../docs/dataset.tsv",
  delim = "\t"
) %>%
  select(
    database,
    type
  )

structureDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary)

organismDictionary <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

## temp path
ott_taxonomy <-
  vroom(file = "../data/external/taxonomySource/organism/taxonomy.tsv")

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
  mutate(curationTypeId_1 = 1)

automaticallyValidated <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationTypeId_2 = 2)

manuallyRemoved <-
  vroom_read_safe(path = "../data/validation/manuallyRemoved.tsv.gz") %>%
  select(colnames(inhouseDbMinimal)) %>%
  mutate(curationTypeId_3 = 3)

## fixed for the moment, will have to change
curation_type <-
  data.frame(
    id = c(1, 2, 3, 4),
    name = c(
      "manual_validation",
      "automatic_validation",
      "manual_removal",
      "automatic_removal"
    )
  )

if (db_mode == "normal") {
  curation_type_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    *
    FROM curation_type"
  )

  curation_type <- curation_type %>%
    anti_join(., curation_type_old) %>%
    distinct(name)

  dbWriteTable(
    conn = db,
    name = "curation_type",
    value = curation_type,
    row.names = FALSE,
    append = TRUE
  )

  curation_type <- dbGetQuery(
    conn = db,
    statement = "SELECT
    *
    FROM curation_type"
  )
}


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

if (db_mode == "fromScratch") {
  dbList <- lapply(pathDataInterimDbDir, vroom_read_safe)

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
}

if (db_mode == "normal") {
  originalTable <-
    vroom_read_safe(path = pathDataInterimTablesOriginalTable)
}

if (db_mode == "fromScratch") {
  database_source <- originalTable %>%
    distinct(database) %>%
    mutate(databaseSourceId = row_number()) %>%
    left_join(., dbTypes) %>%
    group_by(type) %>%
    mutate(databaseTypeId = cur_group_id()) %>%
    ungroup()
}

if (db_mode == "normal") {
  database_type <- dbTypes %>%
    distinct(type)

  database_type_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    database_type.id AS databaseTypeId,
    database_type.name AS type
    FROM database_type"
  )

  database_type_new <- database_type_old %>%
    full_join(., database_type) %>%
    distinct()

  database_type <- database_type_new %>%
    anti_join(., database_type_old) %>%
    select(name = type)

  dbWriteTable(
    conn = db,
    name = "database_type",
    value = database_type,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "normal") {
  database_source <- originalTable %>%
    distinct(database) %>%
    left_join(., dbTypes)

  database_type <- dbGetQuery(
    conn = db,
    statement = "SELECT
    database_type.id AS databaseTypeId,
    database_type.name AS type
    FROM database_type"
  )

  database_source <- database_source %>%
    left_join(., database_type) %>%
    select(
      name = database,
      databaseTypeId
    )

  database_source_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    name
    FROM database_source"
  )

  database_source <- database_source %>%
    anti_join(., database_source_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "database_source",
    value = database_source,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  organism_source <- originalTable %>%
    distinct(organismType, organismValue) %>%
    mutate(organismSourceId = row_number()) %>%
    group_by(organismType) %>%
    mutate(organismTypeId = cur_group_id()) %>%
    ungroup()
}

if (db_mode == "normal") {
  organism_type <- originalTable %>%
    distinct(organismType)

  organism_type_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    organism_type.id AS organismTypeId,
    organism_type.name AS organismType
    FROM organism_type"
  )

  organism_type_new <- organism_type_old %>%
    full_join(., organism_type) %>%
    distinct()

  organism_type <- organism_type_new %>%
    anti_join(., organism_type_old) %>%
    select(name = organismType) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_type",
    value = organism_type,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "normal") {
  organism_type <- dbGetQuery(
    conn = db,
    statement = "SELECT
    organism_type.id AS organismTypeId,
    organism_type.name AS organismType
    FROM organism_type"
  )

  organism_source_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM organism_source"
  )

  organism_source <- originalTable %>%
    distinct(organismType, organismValue) %>%
    left_join(., organism_type) %>%
    select(
      value = organismValue,
      organismTypeId
    ) %>%
    anti_join(., organism_source_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_source",
    value = organism_source,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  structure_source <- originalTable %>%
    distinct(structureType, structureValue) %>%
    mutate(structureSourceId = row_number()) %>%
    group_by(structureType) %>%
    mutate(structureTypeId = cur_group_id()) %>%
    ungroup()
}

if (db_mode == "normal") {
  structure_type <- originalTable %>%
    distinct(structureType)

  structure_type_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    structure_type.id AS structureTypeId,
    structure_type.name AS structureType
    FROM structure_type"
  )

  structure_type_new <- structure_type_old %>%
    full_join(., structure_type) %>%
    distinct()

  structure_type <- structure_type_new %>%
    anti_join(., structure_type_old) %>%
    select(name = structureType) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "structure_type",
    value = structure_type,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "normal") {
  structure_type <- dbGetQuery(
    conn = db,
    statement = "SELECT
    structure_type.id AS structureTypeId,
    structure_type.name AS structureType
    FROM structure_type"
  )

  structure_source_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM structure_source"
  )

  structure_source <- originalTable %>%
    distinct(structureType, structureValue) %>%
    left_join(., structure_type) %>%
    select(
      value = structureValue,
      structureTypeId
    ) %>%
    anti_join(., structure_source_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "structure_source",
    value = structure_source,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  reference_source <- originalTable %>%
    distinct(referenceType, referenceValue) %>%
    mutate(referenceSourceId = row_number()) %>%
    group_by(referenceType) %>%
    mutate(referenceTypeId = cur_group_id()) %>%
    ungroup()
}

if (db_mode == "normal") {
  reference_type <- originalTable %>%
    distinct(referenceType)

  reference_type_old <- dbGetQuery(
    conn = db,
    statement = "SELECT
    reference_type.id AS referenceTypeId,
    reference_type.name AS referenceType
    FROM reference_type"
  )

  reference_type_new <- reference_type_old %>%
    full_join(., reference_type) %>%
    distinct()

  reference_type <- reference_type_new %>%
    anti_join(., reference_type_old) %>%
    select(name = referenceType) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "reference_type",
    value = reference_type,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "normal") {
  reference_type <- dbGetQuery(
    conn = db,
    statement = "SELECT
    reference_type.id AS referenceTypeId,
    reference_type.name AS referenceType
    FROM reference_type"
  )

  reference_source_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM reference_source"
  )

  reference_source <- originalTable %>%
    distinct(referenceType, referenceValue) %>%
    left_join(., reference_type) %>%
    select(
      value = referenceValue,
      referenceTypeId
    ) %>%
    anti_join(., reference_source_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "reference_source",
    value = reference_source,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  database_type <- database_source %>%
    select(
      id = databaseTypeId,
      name = type
    ) %>%
    distinct() %>%
    arrange(id)
}

if (db_mode == "fromScratch") {
  organism_type <- organism_source %>%
    select(
      id = organismTypeId,
      name = organismType
    ) %>%
    distinct() %>%
    arrange(id)
}

if (db_mode == "fromScratch") {
  structure_type <- structure_source %>%
    select(
      id = structureTypeId,
      name = structureType
    ) %>%
    distinct() %>%
    arrange(id)
}

if (db_mode == "fromScratch") {
  reference_type <- reference_source %>%
    select(
      id = referenceTypeId,
      name = referenceType
    ) %>%
    distinct() %>%
    arrange(id)
}

if (db_mode == "fromScratch") {
  data_source <- originalTable %>%
    left_join(., database_source %>% select(-databaseTypeId)) %>%
    left_join(., organism_source %>% select(-organismTypeId)) %>%
    left_join(., structure_source %>% select(-structureTypeId)) %>%
    left_join(., reference_source %>% select(-referenceTypeId)) %>%
    mutate(id = row_number())
}

if (db_mode == "normal") {
  database_source <- dbGetQuery(
    conn = db,
    statement = "
                                SELECT
                                *
                                FROM database_source"
  )

  organism_source <- dbGetQuery(
    conn = db,
    statement = "
                                SELECT
                                *
                                FROM organism_source"
  )

  structure_source <- dbGetQuery(
    conn = db,
    statement = "
                                 SELECT
                                 *
                                 FROM structure_source"
  )

  reference_source <- dbGetQuery(
    conn = db,
    statement = "SELECT
                      *
                      FROM reference_source"
  )

  data_source <- originalTable %>%
    left_join(
      .,
      database_source %>% select(
        database = name,
        databaseSourceId = id,
        everything()
      )
    ) %>%
    left_join(
      .,
      organism_source %>% select(
        organismValue = value,
        organismSourceId = id,
        everything()
      )
    ) %>%
    left_join(
      .,
      structure_source %>% select(
        structureValue = value,
        structureSourceId = id,
        everything()
      )
    ) %>%
    left_join(
      .,
      reference_source %>% select(
        referenceValue = value,
        referenceSourceId = id,
        everything()
      )
    ) %>%
    distinct(
      databaseSourceId,
      organismSourceId,
      structureSourceId,
      referenceSourceId
    )

  data_source_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM data_source"
  )

  data_source <- data_source %>%
    anti_join(., data_source_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "data_source",
    value = data_source,
    row.names = FALSE,
    append = TRUE
  )

  data_source <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM data_source"
  )
}

if (db_mode == "fromScratch") {
  database_source <- database_source %>%
    select(
      id = databaseSourceId,
      name = database,
      databaseTypeId
    )

  organism_source <- organism_source %>%
    select(
      id = organismSourceId,
      value = organismValue,
      organismTypeId
    )

  structure_source <- structure_source %>%
    select(
      id = structureSourceId,
      value = structureValue,
      structureTypeId
    )

  reference_source <- reference_source %>%
    select(
      id = referenceSourceId,
      value = referenceValue,
      referenceTypeId
    )
}

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

if (db_mode == "normal") {
  organism_cleaned_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    name
    FROM organism_cleaned"
  )

  organism_cleaned <- organismOld %>%
    select(name = organismCleaned) %>%
    anti_join(., organism_cleaned_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_cleaned",
    value = organism_cleaned,
    row.names = FALSE,
    append = TRUE
  )

  organism_cleaned <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    id,
    name
    FROM organism_cleaned"
  )
}

if (db_mode == "fromScratch") {
  organism_cleaned <- organismOld %>%
    distinct(organismCleaned) %>%
    mutate(id = row_number()) %>%
    select(id,
      name = organismCleaned
    )
}

organism_synonym <- organism_detected %>%
  left_join(.,
    organism_cleaned,
    by = c("organismCleaned" = "name")
  ) %>%
  select(
    id = id.x,
    name = organismDetected,
    organismCleanedId = id.y
  )

if (db_mode == "normal") {
  organism_synonym_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    name,
    organismCleanedId
    FROM organism_synonym"
  )

  organism_synonym <- organism_synonym %>%
    select(name, organismCleanedId) %>%
    anti_join(., organism_synonym_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_synonym",
    value = organism_synonym,
    row.names = FALSE,
    append = TRUE
  )

  organism_database_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    name
    FROM organism_database"
  )

  organism_database <- organismOld %>%
    distinct(organismCleaned_dbTaxo) %>%
    select(name = organismCleaned_dbTaxo) %>%
    anti_join(., organism_database_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_database",
    value = organism_database,
    row.names = FALSE,
    append = TRUE
  )

  organism_database <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM organism_database"
  )
}

if (db_mode == "fromScratch") {
  organism_database <- organismOld %>%
    distinct(organismCleaned_dbTaxo) %>%
    group_by(organismCleaned_dbTaxo) %>%
    mutate(id = cur_group_id()) %>%
    ungroup() %>%
    select(id,
      name = organismCleaned_dbTaxo
    )
}

if (db_mode == "normal") {
  organism_information_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
      organismCleanedId,
      organismDatabaseId,
      taxonId,
      ranks,
      taxonomy,
      rank
    FROM organism_information"
  )

  organism_information <- organismOld %>%
    left_join(.,
      organism_cleaned,
      by = c("organismCleaned" = "name")
    ) %>%
    select(
      organismCleanedId = id,
      everything()
    ) %>%
    left_join(.,
      organism_database,
      by = c("organismCleaned_dbTaxo" = "name")
    ) %>%
    select(
      organismDatabaseId = id,
      everything()
    ) %>%
    distinct(organismCleanedId,
      organismDatabaseId,
      .keep_all = TRUE
    ) %>%
    anti_join(organism_information_old) %>%
    select(
      organismCleanedId,
      organismDatabaseId,
      taxonId = organismCleaned_id,
      ranks = organismCleaned_dbTaxoTaxonRanks,
      taxonomy = organismCleaned_dbTaxoTaxonomy,
      rank = organismCleaned_rank
    ) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "organism_information",
    value = organism_information,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  organism_information <- organismOld %>%
    left_join(.,
      organism_cleaned,
      by = c("organismCleaned" = "name")
    ) %>%
    select(
      organismCleanedId = id,
      everything()
    ) %>%
    left_join(.,
      organism_database,
      by = c("organismCleaned_dbTaxo" = "name")
    ) %>%
    select(
      organismDatabaseId = id,
      everything()
    ) %>%
    distinct(organismCleanedId,
      organismDatabaseId,
      .keep_all = TRUE
    ) %>%
    mutate(id = row_number()) %>%
    select(
      id,
      organismCleanedId,
      organismDatabaseId,
      taxonId = organismCleaned_id,
      ranks = organismCleaned_dbTaxoTaxonRanks,
      taxonomy = organismCleaned_dbTaxoTaxonomy,
      rank = organismCleaned_rank
    )
}

if (db_mode == "normal") {
  reference_cleaned_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM reference_cleaned"
  )

  reference_cleaned <- referenceOrganismDictionary %>%
    filter(!is.na(referenceCleanedTitle)) %>%
    distinct(
      referenceCleanedDoi,
      referenceCleanedPmcid,
      referenceCleanedPmid,
      referenceCleanedTitle
    ) %>%
    left_join(., referenceOrganismDictionary) %>%
    select(
      doi = referenceCleanedDoi,
      pmcid = referenceCleanedPmcid,
      pmid = referenceCleanedPmid,
      title = referenceCleanedTitle
    ) %>%
    anti_join(., reference_cleaned_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "reference_cleaned",
    value = reference_cleaned,
    row.names = FALSE,
    append = TRUE
  )

  reference_cleaned <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM reference_cleaned"
  )
}

if (db_mode == "fromScratch") {
  reference_cleaned <- referenceOrganismDictionary %>%
    filter(!is.na(referenceCleanedTitle)) %>%
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
}

if (db_mode == "normal") {
  structure_cleaned_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM structure_cleaned"
  )

  structure_cleaned <- structureOld %>%
    distinct(
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
    mutate(
      stereocentersTotal = as.integer(structureCleaned_stereocenters_total),
      stereocentersUnspecified = as.integer(structureCleaned_stereocenters_unspecified),
      exactMass = as.numeric(structureCleaned_exactMass),
      xlogp = as.numeric(structureCleaned_xlogp)
    ) %>%
    select(
      traditionalName = structureCleaned_nameTraditional,
      iupacName = structureCleaned_nameIupac,
      inchikey = structureCleanedInchikey,
      inchikey2D = structureCleaned_inchikey2D,
      inchi = structureCleanedInchi,
      inchi2D = structureCleaned_inchi2D,
      smiles = structureCleanedSmiles,
      smiles2D = structureCleaned_smiles2D,
      stereocentersTotal,
      stereocentersUnspecified,
      molecularFormula = structureCleaned_molecularFormula,
      exactMass,
      xlogp
    ) %>%
    anti_join(., structure_cleaned_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "structure_cleaned",
    value = structure_cleaned,
    row.names = FALSE,
    append = TRUE
  )

  structure_cleaned <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM structure_cleaned"
  )
}

if (db_mode == "fromScratch") {
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
}

if (db_mode == "fromScratch") {
  inhouseDbMinimal_complemented <- inhouseDbMinimal %>%
    mutate(curationTypeId_4 = 4) %>%
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
      curationTypeId = value
    )
}

if (db_mode == "normal") {
  inhouseDbMinimal_complemented <-
    semi_join(inhouseDbMinimal, originalTable) %>%
    mutate(curationTypeId_4 = 4) %>%
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
      curationTypeId = value
    )
}

data_cleaned_temp <- inhouseDbMinimal_complemented %>%
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
  left_join(.,
    reference_cleaned %>% distinct(id, title),
    by = c("referenceCleanedTitle" = "title")
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
    curationTypeId
  )

if (db_mode == "normal") {
  data_cleaned <- data_cleaned_temp %>%
    distinct(
      organismCleanedId,
      structureCleanedId,
      referenceCleanedId,
      curationTypeId
    )

  data_cleaned_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM data_cleaned"
  )

  data_cleaned <- data_cleaned %>%
    anti_join(., data_cleaned_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "data_cleaned",
    value = data_cleaned,
    row.names = FALSE,
    append = TRUE
  )

  data_cleaned <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM data_cleaned"
  )
}

if (db_mode == "normal") {
  data_source__data_cleaned_old <- dbGetQuery(
    conn = db,
    statement = "
    SELECT
    *
    FROM data_source__data_cleaned"
  )

  data_source <- dbGetQuery(
    conn = db,
    statement = sqlFromFile("schema_db/0001_extract_data_source.sql")
  )

  data_source__data_cleaned <- data_cleaned_temp %>%
    left_join(
      .,
      data_source %>%
        select(
          dataSourceId = id,
          everything()
        )
    ) %>%
    left_join(., data_cleaned %>%
      select(
        dataCleanedId = id,
        everything()
      )) %>%
    filter(!is.na(dataCleanedId)) %>%
    distinct(dataSourceId, dataCleanedId)

  data_source__data_cleaned <- data_source__data_cleaned %>%
    anti_join(., data_source__data_cleaned_old) %>%
    distinct()

  dbWriteTable(
    conn = db,
    name = "data_source__data_cleaned",
    value = data_source__data_cleaned,
    row.names = FALSE,
    append = TRUE
  )
}

if (db_mode == "fromScratch") {
  data_cleaned <- data_cleaned_temp %>%
    distinct(
      organismCleanedId,
      structureCleanedId,
      referenceCleanedId,
      curationTypeId
    ) %>%
    mutate(id = row_number()) %>%
    left_join(., data_cleaned_temp)

  data_source__data_cleaned <- data_cleaned_temp %>%
    left_join(
      .,
      data_source %>%
        select(
          dataSourceId = id,
          everything(),
          -type
        )
    ) %>%
    left_join(., data_cleaned %>%
      select(
        dataCleanedId = id,
        everything()
      )) %>%
    filter(!is.na(dataCleanedId)) %>%
    distinct(dataSourceId, dataCleanedId)

  data_source__data_cleaned <- data_source__data_cleaned %>%
    mutate(id = row_number()) %>%
    select(
      id,
      dataSourceId,
      dataCleanedId
    )

  data_cleaned <- data_cleaned %>%
    distinct(
      id,
      structureCleanedId,
      organismCleanedId,
      referenceCleanedId,
      curationTypeId
    ) %>%
    select(
      id,
      structureCleanedId,
      organismCleanedId,
      referenceCleanedId,
      curationTypeId
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
}

if (db_mode == "fromScratch") {
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

  cat("checking if tables match \n")
  identical(
    dbListFields(db, "curation_type"),
    colnames(curation_type)
  )
  identical(
    dbListFields(db, "data_cleaned"),
    colnames(data_cleaned)
  )
  identical(
    dbListFields(db, "data_source"),
    colnames(data_source)
  )
  identical(
    dbListFields(db, "data_source__data_cleaned"),
    colnames(data_source__data_cleaned)
  )
  identical(
    dbListFields(db, "database_source"),
    colnames(database_source)
  )
  identical(
    dbListFields(db, "database_type"),
    colnames(database_type)
  )
  identical(
    dbListFields(db, "organism_cleaned"),
    colnames(organism_cleaned)
  )
  identical(
    dbListFields(db, "organism_database"),
    colnames(organism_database)
  )
  identical(
    dbListFields(db, "organism_information"),
    colnames(organism_information)
  )
  identical(
    dbListFields(db, "organism_source"),
    colnames(organism_source)
  )
  identical(
    dbListFields(db, "organism_synonym"),
    colnames(organism_synonym)
  )
  identical(
    dbListFields(db, "organism_type"),
    colnames(organism_type)
  )
  identical(
    dbListFields(db, "ott_taxonomy"),
    colnames(ott_taxonomy)
  )
  identical(
    dbListFields(db, "reference_cleaned"),
    colnames(reference_cleaned)
  )
  # identical(dbListFields(db, "reference_database"),
  #           colnames(reference_database)) ## TO DO
  # identica(dbListFields(db, "reference_information"),
  #          colnames(reference_information)) ## TO DO
  identical(
    dbListFields(db, "reference_source"),
    colnames(reference_source)
  )
  identical(
    dbListFields(db, "reference_type"),
    colnames(reference_type)
  )
  identical(
    dbListFields(db, "structure_cleaned"),
    colnames(structure_cleaned)
  )
  # identical(
  #   dbListFields(db, "structure_cleaned__structure_information"),
  #   colnames(structure_cleaned__structure_information)
  # ) ## TO DO
  # identical(dbListFields(db, "structure_information"),
  #           colnames(structure_information)) ## TO DO
  identical(
    dbListFields(db, "structure_source"),
    colnames(structure_source)
  )
  identical(
    dbListFields(db, "structure_type"),
    colnames(structure_type)
  )

  dbWriteTable(
    conn = db,
    name = "curation_type",
    value = curation_type,
    row.names = FALSE,
    append = TRUE
  )

  dbWriteTable(
    conn = db,
    name = "data_cleaned",
    value = data_cleaned,
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
    name = "data_source__data_cleaned",
    value = data_source__data_cleaned,
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
    name = "database_type",
    value = database_type,
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
    name = "organism_database",
    value = organism_database,
    row.names = FALSE,
    append = TRUE
  )

  dbWriteTable(
    conn = db,
    name = "organism_information",
    value = organism_information,
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
    name = "organism_synonym",
    value = organism_synonym,
    row.names = FALSE,
    append = TRUE
  )

  dbWriteTable(
    conn = db,
    name = "organism_type",
    value = organism_type,
    row.names = FALSE,
    append = TRUE
  )

  dbWriteTable(
    conn = db,
    name = "ott_taxonomy",
    value = ott_taxonomy,
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
  #   name = "reference_database",
  #   value = reference_database,
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
    name = "reference_type",
    value = reference_type,
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

  # dbWriteTable(
  #   conn = db,
  #   name = "structure_cleaned__structure_information",
  #   value = structure_cleaned__structure_information,
  #   row.names = FALSE,
  #   append = TRUE
  # )

  # dbWriteTable(
  #   conn = db,
  #   name = "structure_information",
  #   value = structure_information,
  #   row.names = FALSE,
  #   append = TRUE
  # )

  dbWriteTable(
    conn = db,
    name = "structure_source",
    value = structure_source,
    row.names = FALSE,
    append = TRUE
  )

  dbWriteTable(
    conn = db,
    name = "structure_type",
    value = structure_type,
    row.names = FALSE,
    append = TRUE
  )
}

dbDisconnect(db)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
