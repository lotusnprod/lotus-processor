cat(
  "This script aligns all DBs with previous results and outputs files \n",
  "containing new entries for editing \n"
)

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/split_data_table_quote.R")
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(DBI)
library(RSQLite)
library(tidyverse)

cat("... files ... \n")
cat("... DBs \n")

if (mode != "test") {
  ## work in progress for SQLite integration
  # cat("... connecting to SQLite database ...")
  # drv <- SQLite()
  # 
  # db <- dbConnect(
  #   drv = drv,
  #   dbname = lotusDB
  # )
  # 
  # cat("... listing remote objects")
  # dbListObjects(db)
  # 
  # test <- RSQLite::dbGetQuery(
  #   conn = db,
  #   statement = sqlFromFile("schema_db/0001_extract_data_source.sql"))
  
  cat("... list of source databases")
  dbList <- lapply(pathDataInterimDbDir, vroom_read_safe)

  cat("... dictionaries ... \n")
  if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
    cat("... structures \n")
    structureDictionary <-
      vroom_read_safe(path = pathDataInterimDictionariesStructureDictionary) %>%
      tibble()
  }

  if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
    cat("... previously unsucessfully querried structures \n")
    structureAntiDictionary <-
      vroom_read_safe(path = pathDataInterimDictionariesStructureAntiDictionary) %>%
      tibble()
  }

  if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
    cat("... organisms \n")
    organismDictionary <-
      vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary) %>%
      tibble()
  }

  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    cat("... references \n")
    referenceDictionary <-
      vroom_read_safe(path = pathDataInterimDictionariesReferenceDictionary) %>%
      data.table()
  }

  cat("renaming and selecting columns \n")
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

  if (mode == "min") {
    cat("sampling rows for min mode \n")
    set.seed(
      seed = 42,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion"
    )
    dbTable <- dbTable %>%
      sample_n(size = 1000)
  }
}

if (mode == "test") {
  dbTable <- vroom_read_safe(path = pathTests) %>%
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
    ) %>%
    tibble()
}

cat("... full original table \n")
originalTable <- dbTable %>%
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

cat("keeping entries not previously curated only ... \n")
cat("... inchi table \n")
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

cat("... smiles table \n")
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

cat("... chemical names table \n")
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

cat("... structures table \n")
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

cat("... organisms tables ... \n")
cat("... clean \n")
organismTable_clean <- originalTable %>%
  filter(organismType == "clean") %>%
  filter(!is.na(organismValue)) %>%
  distinct(organismValue)

cat("... dirty \n")
organismTable_dirty <- originalTable %>%
  filter(organismType == "dirty") %>%
  filter(!is.na(organismValue)) %>%
  distinct(organismValue)


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

cat("... structures table \n")
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

cat("... references tables ... \n")
cat("... DOI table \n")
referenceTable_doi <- dbTable %>%
  filter(!is.na(referenceOriginal_doi)) %>%
  distinct(referenceOriginal_doi) %>%
  select(referenceOriginal = referenceOriginal_doi)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_doi <- anti_join(
      x = referenceTable_doi,
      y = referenceDictionary
    )
  }
}

referenceTable_doi <- referenceTable_doi %>%
  select(referenceOriginal_doi = referenceOriginal)

if (nrow(referenceTable_doi) == 0) {
  referenceTable_doi[1, ] <- NA
}

cat("... PMID table \n")
referenceTable_pubmed <- dbTable %>%
  filter(!is.na(referenceOriginal_pubmed)) %>%
  distinct(referenceOriginal_pubmed) %>%
  select(referenceOriginal = referenceOriginal_pubmed)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_pubmed <- anti_join(
      x = referenceTable_pubmed,
      y = referenceDictionary
    )
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

cat("... reference title table \n")
referenceTable_title <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(!is.na(referenceOriginal_title)) %>%
  distinct(referenceOriginal_title) %>%
  select(referenceOriginal = referenceOriginal_title) %>%
  data.table()

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_title <- anti_join(
      x = referenceTable_title,
      y = referenceDictionary
    )
  }
}

referenceTable_title <- referenceTable_title %>%
  select(referenceOriginal_title = referenceOriginal)

if (nrow(referenceTable_title) == 0) {
  referenceTable_title[1, ] <- NA
}

cat(".. reference publishing details table \n")
referenceTable_publishingDetails <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(!is.na(referenceOriginal_publishingDetails)) %>%
  distinct(referenceOriginal_publishingDetails) %>%
  select(referenceOriginal = referenceOriginal_publishingDetails)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_publishingDetails <-
      anti_join(
        x = referenceTable_publishingDetails,
        y = referenceDictionary
      )
  }
}

referenceTable_publishingDetails <-
  referenceTable_publishingDetails %>%
  select(referenceOriginal_publishingDetails = referenceOriginal)

if (nrow(referenceTable_publishingDetails) == 0) {
  referenceTable_publishingDetails <-
    rbind(referenceTable_publishingDetails, list(NA))
}

cat("... reference split table \n")
referenceTable_split <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(!is.na(referenceOriginal_split)) %>%
  distinct(referenceOriginal_split) %>%
  select(referenceOriginal = referenceOriginal_split)

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_split <-
      anti_join(
        x = referenceTable_split,
        y = referenceDictionary
      )
  }
}

referenceTable_split <- referenceTable_split %>%
  select(referenceOriginal_split = referenceOriginal)

if (nrow(referenceTable_split) == 0) {
  referenceTable_split <-
    rbind(referenceTable_split, list(NA))
}

cat("... original references table \n")
referenceTable_original <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(is.na(referenceOriginal_split)) %>%
  filter(!is.na(referenceOriginal_original)) %>%
  distinct(referenceOriginal_original) %>%
  select(referenceOriginal = referenceOriginal_original) %>%
  data.table()

if (mode != "test") {
  if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
    referenceTable_original <-
      anti_join(
        x = referenceTable_original,
        y = referenceDictionary
      )
  }
}

referenceTable_original <- referenceTable_original %>%
  select(referenceOriginal_original = referenceOriginal)

if (nrow(referenceTable_original) == 0) {
  referenceTable_original <-
    rbind(referenceTable_original, list(NA))
}

cat("... full references table \n")
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


cat("ensuring directories exist ... \n")
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
    dir.create(pathDataInterimTablesOriginalOrganism,
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
    dir.create(pathDataInterimTablesOriginalOrganism,
      showWarnings = FALSE
    )
)

##### structure
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalStructure),
  yes = dir.create(pathDataInterimTablesOriginalStructure),
  no = paste(pathDataInterimTablesOriginalStructure, "exists")
)

cat("exporting ... \n")
if (nrow(organismTable_clean) != 0) {
  vroom_write(
    x = organismTable_clean,
    path = gzfile(
      description = pathDataInterimTablesOriginalOrganismFile,
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
}

cat(pathDataInterimTablesOriginal, "\n")
vroom_write_safe(
  x = organismTable_full,
  path = pathDataInterimTablesOriginalOrganismFull
)

cat(pathDataInterimTablesOriginalReferenceDoi, "\n")
vroom_write_safe(
  x = referenceTable_doi,
  path = pathDataInterimTablesOriginalReferenceDoi
)

cat(pathDataInterimTablesOriginalReferencePubmed, "\n")
vroom_write_safe(
  x = referenceTable_pubmed,
  path = pathDataInterimTablesOriginalReferencePubmed
)

cat(pathDataInterimTablesOriginalReferenceTitleFolder, "\n")
split_data_table_quote(
  x = referenceTable_title,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceTitleFolder
)

cat(
  pathDataInterimTablesOriginalReferencePublishingDetails,
  "\n"
)
vroom_write_safe(
  x = referenceTable_publishingDetails,
  path = pathDataInterimTablesOriginalReferencePublishingDetails
)

cat(pathDataInterimTablesOriginalReferenceSplit, "\n")
vroom_write_safe(
  x = referenceTable_split,
  path = pathDataInterimTablesOriginalReferenceSplit
)

cat(pathDataInterimTablesOriginalReferenceOriginalFolder, "\n")
split_data_table_quote(
  x = referenceTable_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

cat(pathDataInterimTablesOriginalReferenceFull, "\n")
vroom_write_safe(
  x = referenceTable_full,
  path = pathDataInterimTablesOriginalReferenceFull
)

cat(pathDataInterimTablesOriginalStructureInchi, "\n")
vroom_write_safe(
  x = structureTable_inchi,
  path = pathDataInterimTablesOriginalStructureInchi
)

cat(pathDataInterimTablesOriginalStructureNominal, "\n")
vroom_write_safe(
  x = structureTable_nominal,
  path = pathDataInterimTablesOriginalStructureNominal
)

cat(pathDataInterimTablesOriginalStructureSmiles, "\n")
vroom_write_safe(
  x = structureTable_smiles,
  path = pathDataInterimTablesOriginalStructureSmiles
)

cat(pathDataInterimTablesOriginalStructureFull, "\n")
vroom_write_safe(
  x = structureTable_full,
  path = pathDataInterimTablesOriginalStructureFull
)

cat(pathDataInterimTablesOriginalTable, "\n")
vroom_write_safe(
  x = originalTable,
  path = pathDataInterimTablesOriginalTable
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
