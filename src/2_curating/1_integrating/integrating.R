cat(
  "This script aligns all DBs with previous results and outputs files \n",
  "containing new entries for editing \n"
)

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions/helpers.R")

cat("loading ... \n")
cat("... libraries \n")
library(data.table)
library(tidyverse)

cat("... files ... \n")
cat("... DBs \n")
dbList <- lapply(pathDataInterimDbDir, db_loader)

cat("... dictionaries ... \n")
cat("... structures \n")
if (file.exists(pathDataInterimDictionariesStructureDictionary)) {
  structureDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesStructureDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble()
}

cat("... previously unsucessfully querried structures \n")
if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
  structureAntiDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesStructureAntiDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble()
}

cat("... organisms \n")
if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  organismDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesOrganismDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble()
}

cat("... references \n")
if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesReferenceDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    data.table()
}

cat("renaming and selecting columns \n")
dbTable <- rbindlist(l = dbList, fill = TRUE) %>%
  select(
    database,
    organismOriginal = biologicalsource,
    structureOriginal_inchi = inchi,
    structureOriginal_nominal = name,
    structureOriginal_smiles = smiles,
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
}
if (mode == "min") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
}
if (mode == "min") {
  dbTable <- dbTable %>%
    sample_n(size = 2000)
}

cat("removing unfriendly characters \n")
dbTable[] <- lapply(dbTable, function(x) {
  gsub("\r\n", " ", x)
})
dbTable[] <- lapply(dbTable, function(x) {
  gsub("\r", " ", x)
})
dbTable[] <- lapply(dbTable, function(x) {
  gsub("\n", " ", x)
})
dbTable[] <- lapply(dbTable, function(x) {
  gsub("\t", " ", x)
})

cat("keeping entries not previously curated only ... \n")
cat("... inchi table \n")
structureTable_inchi <- dbTable %>%
  filter(grepl(
    pattern = "^InChI=.*",
    x = structureOriginal_inchi
  )) %>%
  distinct(structureOriginal_inchi) %>%
  select(structureValue = structureOriginal_inchi)

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

structureTable_inchi <- structureTable_inchi %>%
  select(structureOriginal_inchi = structureValue)

if (nrow(structureTable_inchi) == 0) {
  structureTable_inchi <- rbind(structureTable_inchi, list(NA))
}

cat("... smiles table \n")
structureTable_smiles <- dbTable %>%
  filter(!grepl(
    pattern = "^InChI=.*",
    x = structureOriginal_inchi
  )) %>%
  filter(!is.na(structureOriginal_smiles)) %>%
  distinct(structureOriginal_smiles) %>%
  select(structureValue = structureOriginal_smiles)

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

structureTable_smiles <- structureTable_smiles %>%
  select(structureOriginal_smiles = structureValue)

if (nrow(structureTable_smiles) == 0) {
  structureTable_smiles <- rbind(structureTable_smiles, list(NA))
}

cat("... chemical names table \n")
structureTable_nominal <- dbTable %>%
  filter(!grepl(
    pattern = "^InChI=.*",
    x = structureOriginal_inchi
  )) %>%
  filter(is.na(structureOriginal_smiles)) %>%
  filter(!is.na(structureOriginal_nominal)) %>%
  distinct(structureOriginal_nominal) %>%
  select(structureValue = structureOriginal_nominal)

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

structureTable_nominal <- structureTable_nominal %>%
  select(structureOriginal_nominal = structureValue)

if (nrow(structureTable_nominal) == 0) {
  structureTable_nominal <- rbind(structureTable_nominal, list(NA))
}

cat("... structures table \n")
if (file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
  structureTable_full <- dbTable %>%
    filter(!structureOriginal_inchi %in% structureAntiDictionary$structureValue) %>%
    filter(!structureOriginal_smiles %in% structureAntiDictionary$structureValue) %>%
    filter(!structureOriginal_nominal %in% structureAntiDictionary$structureValue) %>%
    filter(!structureOriginal_inchi %in% structureDictionary$structureValue) %>%
    filter(!structureOriginal_smiles %in% structureDictionary$structureValue) %>%
    filter(!structureOriginal_nominal %in% structureDictionary$structureValue)
}

if (!file.exists(pathDataInterimDictionariesStructureAntiDictionary)) {
  structureTable_full <- dbTable
}

structureTable_full <- structureTable_full %>%
  distinct(
    structureOriginal_inchi,
    structureOriginal_smiles,
    structureOriginal_nominal
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 1:ncol(.),
    names_to = c("drop", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(structureType, structureValue)

if (nrow(structureTable_full) == 0) {
  structureTable_full[1, ] <- NA
}

cat("... organisms table \n")
organismTable <- dbTable %>%
  filter(!is.na(organismOriginal)) %>%
  distinct(organismOriginal) %>%
  data.table()

if (file.exists(pathDataInterimDictionariesOrganismDictionary)) {
  organismTable <- anti_join(
    x = organismTable,
    y = organismDictionary
  )
}

cat("... references tables ... \n")
cat("... DOI table \n")
referenceTable_doi <- dbTable %>%
  filter(!is.na(referenceOriginal_doi)) %>%
  distinct(referenceOriginal_doi) %>%
  select(referenceOriginal = referenceOriginal_doi)

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_doi <- anti_join(
    x = referenceTable_doi,
    y = referenceDictionary
  )
}

referenceTable_doi <- referenceTable_doi %>%
  select(referenceOriginal_doi = referenceOriginal)

row.names(referenceTable_doi) <-
  referenceTable_doi$referenceOriginal_doi

cat("... PMID table \n")
referenceTable_pubmed <- dbTable %>%
  filter(!is.na(referenceOriginal_pubmed)) %>%
  distinct(referenceOriginal_pubmed) %>%
  select(referenceOriginal = referenceOriginal_pubmed)

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_pubmed <- anti_join(
    x = referenceTable_pubmed,
    y = referenceDictionary
  )
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

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_title <- anti_join(
    x = referenceTable_title,
    y = referenceDictionary
  )
}

referenceTable_title <- referenceTable_title %>%
  select(referenceOriginal_title = referenceOriginal)

cat(".. reference publishing details table \n")
referenceTable_publishingDetails <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(!is.na(referenceOriginal_publishingDetails)) %>%
  distinct(referenceOriginal_publishingDetails) %>%
  select(referenceOriginal = referenceOriginal_publishingDetails)

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_publishingDetails <-
    anti_join(
      x = referenceTable_publishingDetails,
      y = referenceDictionary
    )
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

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_split <-
    anti_join(
      x = referenceTable_split,
      y = referenceDictionary
    )
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

if (file.exists(pathDataInterimDictionariesReferenceDictionary)) {
  referenceTable_original <-
    anti_join(
      x = referenceTable_original,
      y = referenceDictionary
    )
}

referenceTable_original <- referenceTable_original %>%
  select(referenceOriginal_original = referenceOriginal)

if (nrow(referenceTable_original) == 0) {
  referenceTable_original <-
    rbind(referenceTable_original, list(NA))
}

cat("... full references table \n")
referenceTable_full <- dbTable %>%
  distinct(
    organismOriginal,
    referenceOriginal_authors,
    referenceOriginal_doi,
    referenceOriginal_external,
    referenceOriginal_isbn,
    referenceOriginal_journal,
    referenceOriginal_original,
    referenceOriginal_pubmed,
    referenceOriginal_publishingDetails,
    referenceOriginal_title,
    referenceOriginal_split
  ) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  select(
    organismOriginal,
    referenceType,
    referenceValue
  )

cat("... full original table \n")
originalTable <- dbTable %>%
  select(database, organismOriginal, everything()) %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = c("drop", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = 3:5,
    names_to = c("drop2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(-drop, -drop2)

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
cat(pathDataInterimTablesOriginalOrganism, "\n")
if (nrow(organismTable) != 0) {
  split_data_table(
    x = organismTable,
    no_rows_per_frame = 10000,
    text = "",
    path_to_store = pathDataInterimTablesOriginalOrganism
  )
}

cat(pathDataInterimTablesOriginalReferenceDoi, "\n")
write.table(
  x = referenceTable_doi,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceDoi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalReferencePubmed, "\n")
write.table(
  x = referenceTable_pubmed,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferencePubmed,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = TRUE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
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
write.table(
  x = referenceTable_publishingDetails,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferencePublishingDetails,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalReferenceSplit, "\n")
write.table(
  x = referenceTable_split,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceSplit,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalReferenceOriginalFolder, "\n")
split_data_table_quote(
  x = referenceTable_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

cat(pathDataInterimTablesOriginalReferenceFull, "\n")
write.table(
  x = referenceTable_full,
  file = gzfile(
    description = pathDataInterimTablesOriginalReferenceFull,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalStructureInchi, "\n")
write.table(
  x = structureTable_inchi,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureInchi,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalStructureNominal, "\n")
write.table(
  x = structureTable_nominal,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureNominal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalStructureSmiles, "\n")
write.table(
  x = structureTable_smiles,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureSmiles,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalStructureFull, "\n")
write.table(
  x = structureTable_full,
  file = gzfile(
    description = pathDataInterimTablesOriginalStructureFull,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = TRUE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesOriginalTable, "\n")
write.table(
  x = originalTable,
  file = gzfile(
    description = pathDataInterimTablesOriginalTable,
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