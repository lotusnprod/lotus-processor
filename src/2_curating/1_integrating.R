source("r/log_debug.R")
log_debug(
  "This script aligns all sources with previous results and outputs files \n",
  "containing new entries for editing"
)

start <- Sys.time()
Sys.setlocale("LC_ALL", 'en_US.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... functions")
source("r/split_data_table_quote.R")
source("r/sqlFromFile.R")
source("r/standardizing_original.R")

log_debug("loading ...")
log_debug("... libraries")
library(data.table)
library(DBI)
library(dplyr)
library(readr)
library(tidyr)

log_debug("... files ...")
log_debug("... DBs")

if (mode == "full" | mode == "manual") {
  if (ssot_access == TRUE) {
    library(RPostgreSQL)

    drv <- PostgreSQL()

    log_debug("... connecting to the database")
    # db <- dbConnect(
    #   drv = drv,
    #   dbname = "lotus",
    #   user = "rutza",
    #   host = "localhost",
    # )

    db <- dbConnect(
      drv = drv,
      dbname = dbname,
      user = user,
      host = host,
      port = port,
      password = password
    )

    log_debug("... listing remote objects")
    dbListObjects(db)

    log_debug("... extracting already processed data")
    oldTable <- dbGetQuery(
      conn = db,
      statement = sqlFromFile("queries_db/extract_data_source.sql")
    )
  } else {
    oldTable <- data.frame() %>%
      mutate(
        database = NA,
        organism_value = NA,
        organism_type = NA,
        reference_value = NA,
        reference_type = NA,
        structure_value = NA,
        structure_type = NA
      ) %>%
      mutate_all(as.character)
  }
  log_debug("... list of source databases")
  dbList <- lapply(
    pathDataInterimDbDir,
    read_delim,
    delim = "\t",
    col_types = cols(.default = "c")
  )

  log_debug("... dictionaries ...")
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    log_debug("... structures")
    structureDictionary <-
      read_delim(
        file = pathDataInterimDictionariesStructureDictionary,
        delim = "\t",
        col_types = cols(.default = "c")
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    log_debug("... previously unsucessfully querried structures")
    structureAntiDictionary <-
      read_delim(
        file = pathDataInterimDictionariesStructureAntiDictionary,
        delim = "\t",
        col_types = cols(.default = "c")
      )
  }

  if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
    log_debug("... organisms")
    organismDictionary <-
      read_delim(
        file = pathDataInterimDictionariesOrganismDictionary,
        delim = "\t",
        col_types = cols(.default = "c")
      )
  }

  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    log_debug("... references")
    referenceDictionary <-
      read_delim(
        file = pathDataInterimDictionariesReferenceDictionary,
        delim = "\t",
        col_types = cols(.default = "c")
      ) %>%
      data.table()
  }

  log_debug("renaming and selecting columns")
  dbTable <- rbindlist(l = dbList, fill = TRUE) %>%
    data.frame()

  dbTable[setdiff(
    x = accepted_fields,
    y = names(dbTable)
  )] <- NA

  dbTable <- dbTable %>%
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
    mutate(referenceOriginal_pubmed = as.character(referenceOriginal_pubmed))

  if (mode == "full") {
    log_debug("sampling rows for test mode")
    "%ni%" <- Negate("%in%")
    set.seed(
      seed = 42,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion"
    )
    dbTable_sampled_1 <- dbTable %>%
      filter(database %ni% forbidden_export) %>%
      filter(is.na(referenceOriginal_external)) %>%
      filter(!(
        !is.na(structureOriginal_inchi) &
          !is.na(structureOriginal_nominal)
      ) | !(
        !is.na(structureOriginal_smiles) &
          !is.na(structureOriginal_nominal)
      )) %>% # to avoid too many names (long for CI)
      sample_n(size = 490)
    set.seed(
      seed = 42,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion"
    )
    dbTable_sampled_2 <- dbTable %>%
      filter(database %ni% forbidden_export) %>%
      filter(!is.na(organismOriginal_dirty)) %>%
      filter(is.na(referenceOriginal_external)) %>%
      filter(!(
        !is.na(structureOriginal_inchi) &
          !is.na(structureOriginal_nominal)
      ) | !(
        !is.na(structureOriginal_smiles) &
          !is.na(structureOriginal_nominal)
      )) %>% # to avoid too many names (long for CI)
      sample_n(size = 10)
    dbTable_sampled <- bind_rows(dbTable_sampled_1, dbTable_sampled_2)
    originalTable_sampled <- dbTable_sampled %>%
      select(database, everything()) %>%
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
} else {
  dbTable <- read_delim(
    file = pathTestsFile,
    delim = "\t",
    col_types = cols(.default = "c")
  ) %>%
    pivot_wider(
      names_from = "organismType",
      values_from = "organismValue",
      names_prefix = "organismOriginal_"
    ) %>%
    pivot_wider(
      names_from = "structureType",
      values_from = "structureValue",
      names_prefix = "structureOriginal_"
    ) %>%
    pivot_wider(
      names_from = "referenceType",
      values_from = "referenceValue",
      names_prefix = "referenceOriginal_"
    ) %>%
    unnest() %>%
    mutate(
      referenceOriginal_authors = NA,
      referenceOriginal_external = NA,
      referenceOriginal_isbn = NA,
      referenceOriginal_journal = NA
    ) %>%
    select(
      database,
      organismOriginal_clean,
      organismOriginal_dirty,
      structureOriginal_inchi,
      structureOriginal_nominal,
      structureOriginal_smiles,
      referenceOriginal_authors,
      referenceOriginal_doi,
      referenceOriginal_external,
      referenceOriginal_isbn,
      referenceOriginal_journal,
      referenceOriginal_original,
      referenceOriginal_pubmed,
      referenceOriginal_publishingDetails,
      referenceOriginal_split,
      referenceOriginal_title,
    )
}

log_debug("... full original table")
originalTable <- dbTable %>%
  select(database, everything()) %>%
  pivot_longer(
    cols = 7:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop) %>%
  pivot_longer(
    cols = 4:6,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop2) %>%
  pivot_longer(
    cols = 2:3,
    names_to = c("drop3", "organismType"),
    names_sep = "_",
    values_to = "organismValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop3) %>%
  distinct()

if (mode == "full" | mode == "manual") {
  log_debug("new entries only ...")
  originalTable <- anti_join(
    originalTable,
    oldTable %>%
      select(
        database,
        organismValue = organism_value,
        organismType = organism_type,
        referenceValue = reference_value,
        referenceType = reference_type,
        structureValue = structure_value,
        structureType = structure_type
      )
  )
}

log_debug("ensuring proper encoding ...")
originalTable$organismValue <-
  iconv(
    x = originalTable$organismValue,
    to = "UTF-8",
    sub = "Unicode"
  )

originalTable$referenceValue <-
  iconv(
    x = originalTable$referenceValue,
    to = "UTF-8",
    sub = "Unicode"
  )

originalTable$structureValue <-
  iconv(
    x = originalTable$structureValue,
    to = "UTF-8",
    sub = "Unicode"
  )

log_debug("keeping entries not previously curated only ...")
log_debug("... inchi table")
structureTable_inchi <- originalTable %>%
  filter(structureType == "inchi") %>%
  filter(!is.na(structureValue)) %>%
  distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_inchi <-
      anti_join(
        x = structureTable_inchi,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_inchi <-
      anti_join(
        x = structureTable_inchi,
        y = structureAntiDictionary
      )
  }
}

structureTable_inchi <- structureTable_inchi %>%
  select(structureOriginal_inchi = structureValue)

if (nrow(structureTable_inchi) == 0) {
  structureTable_inchi <- data.frame(structureOriginal_inchi = NA)
}

log_debug("... smiles table")
structureTable_smiles <- originalTable %>%
  filter(structureType == "smiles") %>%
  filter(!is.na(structureValue)) %>%
  distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_smiles <-
      anti_join(
        x = structureTable_smiles,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_smiles <-
      anti_join(
        x = structureTable_smiles,
        y = structureAntiDictionary
      )
  }
}

structureTable_smiles <- structureTable_smiles %>%
  select(structureOriginal_smiles = structureValue)

if (nrow(structureTable_smiles) == 0) {
  structureTable_smiles <- data.frame(structureOriginal_smiles = NA)
}

log_debug("... chemical names table")
structureTable_nominal <- originalTable %>%
  filter(structureType == "nominal") %>%
  filter(!is.na(structureValue)) %>%
  distinct(structureValue)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    structureTable_nominal <-
      anti_join(
        x = structureTable_nominal,
        y = structureDictionary
      )
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    structureTable_nominal <-
      anti_join(
        x = structureTable_nominal,
        y = structureAntiDictionary
      )
  }
}

structureTable_nominal <- structureTable_nominal %>%
  select(structureOriginal_nominal = structureValue)

if (nrow(structureTable_nominal) == 0) {
  structureTable_nominal <- data.frame(structureOriginal_nominal = NA)
}

log_debug("... structures table")
structureTable_full <-
  bind_rows(
    structureTable_inchi %>%
      mutate(structureType = "inchi") %>%
      select(structureType,
        structureValue = structureOriginal_inchi
      ),
    structureTable_smiles %>%
      mutate(structureType = "smiles") %>%
      select(structureType,
        structureValue = structureOriginal_smiles
      ),
    structureTable_nominal %>%
      mutate(structureType = "nominal") %>%
      select(structureType,
        structureValue = structureOriginal_nominal
      )
  ) %>%
  distinct()

if (nrow(structureTable_full) == 0) {
  structureTable_full[1, ] <- NA
}

log_debug("... organisms tables ...")
log_debug("... clean")
organismTable_clean <- originalTable %>%
  filter(organismType == "clean") %>%
  filter(!is.na(organismValue)) %>%
  distinct(organismValue)

if (nrow(organismTable_clean) == 0) {
  organismTable_clean[1, "organismValue"] <- NA
}

log_debug("... dirty")
organismTable_dirty <- originalTable %>%
  filter(organismType == "dirty") %>%
  filter(!is.na(organismValue)) %>%
  distinct(organismValue)

if (nrow(organismTable_dirty) == 0) {
  organismTable_dirty[1, "organismValue"] <- NA
}

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
    organismTable_clean <- anti_join(
      x = organismTable_clean,
      y = organismDictionary
    )

    organismTable_dirty <- anti_join(
      x = organismTable_dirty,
      y = organismDictionary
    )
  }
}

organismTable_clean <- organismTable_clean %>%
  select(organismOriginal_clean = organismValue)

organismTable_dirty <- organismTable_dirty %>%
  select(organismOriginal_dirty = organismValue) %>%
  data.table()

log_debug("... structures table")
organismTable_full <- bind_rows(
  organismTable_clean %>%
    mutate(organismType = "clean") %>%
    select(organismType,
      organismValue = organismOriginal_clean
    ),
  organismTable_dirty %>%
    mutate(organismType = "dirty") %>%
    select(organismType,
      organismValue = organismOriginal_dirty
    )
) %>%
  distinct()

log_debug("... references tables ...")
log_debug("... DOI table")
referenceTable_doi <- dbTable %>%
  filter(!is.na(referenceOriginal_doi)) %>%
  distinct(referenceOriginal_doi) %>%
  select(referenceOriginal = referenceOriginal_doi) %>%
  mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_doi <- anti_join(
      x = referenceTable_doi %>%
        mutate(origin = "doi"),
      y = referenceDictionary
    ) %>%
      select(-origin)
  }
}

referenceTable_doi <- referenceTable_doi %>%
  select(referenceOriginal_doi = referenceOriginal)

if (nrow(referenceTable_doi) == 0) {
  referenceTable_doi[1, ] <- NA
}

log_debug("... PMID table")
referenceTable_pubmed <- dbTable %>%
  filter(!is.na(referenceOriginal_pubmed)) %>%
  distinct(referenceOriginal_pubmed) %>%
  select(referenceOriginal = referenceOriginal_pubmed) %>%
  mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_pubmed <- anti_join(
      x = referenceTable_pubmed %>%
        mutate(origin = "pubmed"),
      y = referenceDictionary
    ) %>%
      select(-origin)
  }
}

referenceTable_pubmed <- referenceTable_pubmed %>%
  select(referenceOriginal_pubmed = referenceOriginal)

if (nrow(referenceTable_pubmed) == 0) {
  referenceTable_pubmed <- rbind(referenceTable_pubmed, list(NA))
}

if (nrow(referenceTable_pubmed) != 1) {
  row.names(referenceTable_pubmed) <-
    referenceTable_pubmed$referenceOriginal_pubmed
}

log_debug("... reference title table")
referenceTable_title <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(!is.na(referenceOriginal_title)) %>%
  distinct(referenceOriginal_title) %>%
  select(referenceOriginal = referenceOriginal_title) %>%
  mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_title <- anti_join(
      x = referenceTable_title %>%
        mutate(origin = "title"),
      y = referenceDictionary
    ) %>%
      select(-origin)
  }
}

referenceTable_title <- referenceTable_title %>%
  select(referenceOriginal_title = referenceOriginal)

if (nrow(referenceTable_title) == 0) {
  referenceTable_title[1, ] <- NA
}

referenceTable_title <- referenceTable_title %>%
  data.table()

log_debug(".. reference publishing details table")
referenceTable_publishingDetails <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(!is.na(referenceOriginal_publishingDetails)) %>%
  distinct(referenceOriginal_publishingDetails) %>%
  select(referenceOriginal = referenceOriginal_publishingDetails) %>%
  mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_publishingDetails <-
      anti_join(
        x = referenceTable_publishingDetails %>%
          mutate(origin = "publishingDetails"),
        y = referenceDictionary
      ) %>%
      select(-origin)
  }
}

referenceTable_publishingDetails <-
  referenceTable_publishingDetails %>%
  select(referenceOriginal_publishingDetails = referenceOriginal)

if (nrow(referenceTable_publishingDetails) == 0) {
  referenceTable_publishingDetails <-
    rbind(referenceTable_publishingDetails, list(NA))
}

referenceTable_publishingDetails <- referenceTable_publishingDetails %>%
  data.table()

log_debug("... reference split table")
referenceTable_split <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(!is.na(referenceOriginal_split)) %>%
  distinct(referenceOriginal_split) %>%
  select(referenceOriginal = referenceOriginal_split) %>%
  mutate_all(as.character)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_split <-
      anti_join(
        x = referenceTable_split %>%
          mutate(origin = "split"),
        y = referenceDictionary
      ) %>%
      select(-origin)
  }
}

referenceTable_split <- referenceTable_split %>%
  select(referenceOriginal_split = referenceOriginal)

if (nrow(referenceTable_split) == 0) {
  referenceTable_split <-
    rbind(referenceTable_split, list(NA))
}

referenceTable_split <- referenceTable_split %>%
  data.table()

log_debug("... original references table")
referenceTable_original <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(is.na(referenceOriginal_split)) %>%
  filter(!is.na(referenceOriginal_original)) %>%
  distinct(referenceOriginal_original) %>%
  select(referenceOriginal = referenceOriginal_original) %>%
  mutate_all(as.character) %>%
  data.table()

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_original <-
      anti_join(
        x = referenceTable_original %>%
          mutate(origin = "original"),
        y = referenceDictionary
      ) %>%
      select(-origin)
  }
}

referenceTable_original <- referenceTable_original %>%
  select(referenceOriginal_original = referenceOriginal)

if (nrow(referenceTable_original) == 0) {
  referenceTable_original <-
    rbind(referenceTable_original, list(NA))
}

log_debug("... full references table")
referenceTable_full <- originalTable %>%
  select(
    organismType,
    organismValue,
    referenceType,
    referenceValue
  ) %>%
  distinct()

# data_source <- dbGetQuery(
#   conn = db,
#   statement = "SELECT * FROM data_source"
# )

log_debug("ensuring directories exist ...")
ifelse(
  test = !dir.exists(pathDataInterim),
  yes = dir.create(pathDataInterim),
  no = paste(pathDataInterim, "exists")
)
ifelse(
  test = !dir.exists(pathDataInterimTables),
  yes = dir.create(pathDataInterimTables),
  no = paste(pathDataInterimTables, "exists")
)

#### original
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginal),
  yes = dir.create(pathDataInterimTablesOriginal),
  no = paste(pathDataInterimTablesOriginal, "exists")
)

##### organism
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalOrganism),
  yes = dir.create(pathDataInterimTablesOriginalOrganism),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesOriginalOrganism,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesOriginalOrganism,
      showWarnings = FALSE
    )
)

##### reference
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReference),
  yes = dir.create(pathDataInterimTablesOriginalReference),
  no = paste(pathDataInterimTablesOriginalReference, "exists")
)

##### reference folders
###### original
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReferenceOriginalFolder),
  yes = dir.create(pathDataInterimTablesOriginalReferenceOriginalFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesOriginalReferenceOriginalFolder,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesOriginalReferenceOriginalFolder,
      showWarnings = FALSE
    )
)

###### publishing details
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReferencePublishingDetailsFolder),
  yes = dir.create(pathDataInterimTablesOriginalReferencePublishingDetailsFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesOriginalReferencePublishingDetailsFolder,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesOriginalReferencePublishingDetailsFolder,
      showWarnings = FALSE
    )
)

###### split
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReferenceSplitFolder),
  yes = dir.create(pathDataInterimTablesOriginalReferenceSplitFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesOriginalReferenceSplitFolder,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesOriginalReferenceSplitFolder,
      showWarnings = FALSE
    )
)

###### title
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReferenceTitleFolder),
  yes = dir.create(pathDataInterimTablesOriginalReferenceTitleFolder),
  no = file.remove(
    list.files(
      path = pathDataInterimTablesOriginalReferenceTitleFolder,
      full.names = TRUE
    )
  ) &
    dir.create(pathDataInterimTablesOriginalReferenceTitleFolder,
      showWarnings = FALSE
    )
)

##### structure
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalStructure),
  yes = dir.create(pathDataInterimTablesOriginalStructure),
  no = paste(pathDataInterimTablesOriginalStructure, "exists")
)

log_debug("exporting ...")
if (nrow(organismTable_clean) != 0) {
  write_delim(
    x = organismTable_clean,
    delim = "\t",
    file = gzfile(
      description = pathDataInterimTablesOriginalOrganismFile,
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    quote = "none"
  )
  ## because gnverify does not parse quotes
}

log_debug(pathDataInterimTablesOriginal)
write_delim(
  x = organismTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalOrganismFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferenceDoi)
write_delim(
  x = referenceTable_doi,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferenceDoi,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferenceOriginalFolder)
split_data_table_quote(
  x = referenceTable_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

log_debug(pathDataInterimTablesOriginalReferencePubmed)
write_delim(
  x = referenceTable_pubmed,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferencePubmed,
  na = ""
)

log_debug(pathDataInterimTablesOriginalReferencePublishingDetailsFolder)
split_data_table_quote(
  x = referenceTable_publishingDetails,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferencePublishingDetailsFolder
)

log_debug(pathDataInterimTablesOriginalReferenceSplitFolder)
split_data_table_quote(
  x = referenceTable_split,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceSplitFolder
)

log_debug(pathDataInterimTablesOriginalReferenceTitleFolder)
split_data_table_quote(
  x = referenceTable_title,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceTitleFolder
)

log_debug(pathDataInterimTablesOriginalReferenceFull)
write_delim(
  x = referenceTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalReferenceFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureInchi)
write_delim(
  x = structureTable_inchi,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureInchi,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureNominal)
write_delim(
  x = structureTable_nominal,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureNominal,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureSmiles)
write_delim(
  x = structureTable_smiles,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureSmiles,
  na = ""
)

log_debug(pathDataInterimTablesOriginalStructureFull)
write_delim(
  x = structureTable_full,
  delim = "\t",
  file = pathDataInterimTablesOriginalStructureFull,
  na = ""
)

log_debug(pathDataInterimTablesOriginalTable)
write_delim(
  x = originalTable,
  delim = "\t",
  file = pathDataInterimTablesOriginalTable,
  na = ""
)

if (mode == "full") {
  log_debug(pathTestsFile)
  write_delim(
    x = originalTable_sampled,
    delim = "\t",
    file = pathTestsFile,
    na = ""
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
