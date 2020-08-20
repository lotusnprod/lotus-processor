# title: "integrating DBs"

# loading
## paths
source("paths.R")

source("functions/helpers.R")

## libraries
library(data.table)
library(tidyverse)

## files
print(x = "loading DBs")
dbList <- lapply(pathDataInterimDbDir, db_loader)

## dictionaries
print(x = "loading dictionaries")
### structure
if (file.exists(pathDataInterimDictionariesStructureDictionary))
  structureDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesStructureDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

### organism
if (file.exists(pathDataInterimDictionariesOrganismDictionary))
  organismDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesOrganismDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

### reference
if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceDictionary <- read_delim(
    file = gzfile(description = pathDataInterimDictionariesReferenceDictionary),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

# selecting
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
  )

# sampling rows for minimal mode
if (mode == "min")
  set.seed(seed = 42,
           kind = "Mersenne-Twister",
           normal.kind = "Inversion")
if (mode == "min")
  dbTable <- dbTable %>%
  sample_n(size = 2000)

# removing unfriendly characters
dbTable[] <- lapply(dbTable, function(x)
  gsub("\r\n", " ", x))
dbTable[] <- lapply(dbTable, function(x)
  gsub("\r", " ", x))
dbTable[] <- lapply(dbTable, function(x)
  gsub("\n", " ", x))
dbTable[] <- lapply(dbTable, function(x)
  gsub("\t", " ", x))

# sub-setting
print(x = "keeping new data only")

# structures
## with InChI
structureTable_inchi <- dbTable %>%
  filter(grepl(pattern = "^InChI=.*",
               x = structureOriginal_inchi)) %>%
  distinct(structureOriginal_inchi) %>%
  select(structureValue = structureOriginal_inchi)

if (file.exists(pathDataInterimDictionariesStructureDictionary))
  structureTable_inchi <-
  anti_join(x = structureTable_inchi,
            y = structureDictionary)

structureTable_inchi <- structureTable_inchi %>%
  select(structureOriginal_inchi = structureValue)

### without InChI but SMILES
structureTable_smiles <- dbTable %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginal_inchi)) %>%
  filter(!is.na(structureOriginal_smiles)) %>%
  distinct(structureOriginal_smiles) %>%
  select(structureValue = structureOriginal_smiles)

if (file.exists(pathDataInterimDictionariesStructureDictionary))
  structureTable_smiles <-
  anti_join(x = structureTable_smiles,
            y = structureDictionary)

structureTable_smiles <- structureTable_smiles %>%
  select(structureOriginal_smiles = structureValue)

### without InChI nor SMILES but name
structureTable_nominal <- dbTable %>%
  filter(!grepl(pattern = "^InChI=.*",
                x = structureOriginal_inchi)) %>%
  filter(is.na(structureOriginal_smiles)) %>%
  filter(!is.na(structureOriginal_nominal)) %>%
  distinct(structureOriginal_nominal) %>%
  select(structureValue = structureOriginal_nominal)

if (file.exists(pathDataInterimDictionariesStructureDictionary))
  structureTable_nominal <-
  anti_join(x = structureTable_nominal,
            y = structureDictionary)

structureTable_nominal <- structureTable_nominal %>%
  select(structureOriginal_nominal = structureValue)

### full
structureTable_full <- dbTable %>%
  distinct(structureOriginal_inchi,
           structureOriginal_smiles,
           structureOriginal_nominal) %>%
  mutate_all(as.character) %>%
  pivot_longer(
    cols = 1:ncol(.),
    names_to = c("drop", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  select(structureType, structureValue)

## organism
organismTable <- dbTable %>%
  filter(!is.na(organismOriginal)) %>%
  distinct(organismOriginal)

if (file.exists(pathDataInterimDictionariesOrganismDictionary))
  organismTable <- anti_join(x = organismTable,
                             y = organismDictionary)

## reference
### DOI
referenceTable_doi <- dbTable %>%
  filter(!is.na(referenceOriginal_doi)) %>%
  distinct(referenceOriginal_doi) %>%
  select(referenceOriginal = referenceOriginal_doi)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_doi <- anti_join(x = referenceTable_doi,
                                  y = referenceDictionary)

referenceTable_doi <- referenceTable_doi %>%
  select(referenceOriginal_doi = referenceOriginal)


row.names(referenceTable_doi) <-
  referenceTable_doi$referenceOriginal_doi

### pubmed
referenceTable_pubmed <- dbTable %>%
  filter(!is.na(referenceOriginal_pubmed)) %>%
  distinct(referenceOriginal_pubmed) %>%
  select(referenceOriginal = referenceOriginal_pubmed)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_pubmed <- anti_join(x = referenceTable_pubmed,
                                     y = referenceDictionary)

referenceTable_pubmed <- referenceTable_pubmed %>%
  select(referenceOriginal_pubmed = referenceOriginal)

if (nrow(referenceTable_pubmed) == 0)
  referenceTable_pubmed <- rbind(referenceTable_pubmed, list(NA))

if (nrow(referenceTable_pubmed) != 1)
  row.names(referenceTable_pubmed) <-
  referenceTable_pubmed$referenceOriginal_pubmed

### title
referenceTable_title <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(!is.na(referenceOriginal_title)) %>%
  distinct(referenceOriginal_title) %>%
  select(referenceOriginal = referenceOriginal_title)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_title <- anti_join(x = referenceTable_title,
                                    y = referenceDictionary)

referenceTable_title  <- referenceTable_title %>%
  select(referenceOriginal_title = referenceOriginal)

### publishing details
referenceTable_publishingDetails <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(!is.na(referenceOriginal_publishingDetails)) %>%
  distinct(referenceOriginal_publishingDetails) %>%
  select(referenceOriginal = referenceOriginal_publishingDetails)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_publishingDetails <-
  anti_join(x = referenceTable_publishingDetails,
            y = referenceDictionary)

referenceTable_publishingDetails <-
  referenceTable_publishingDetails %>%
  select(referenceOriginal_publishingDetails = referenceOriginal)

if (nrow(referenceTable_publishingDetails) == 0)
  referenceTable_publishingDetails <-
  rbind(referenceTable_publishingDetails, list(NA))

### split
referenceTable_split <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(!is.na(referenceOriginal_split)) %>%
  distinct(referenceOriginal_split) %>%
  select(referenceOriginal = referenceOriginal_split)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_split <-
  anti_join(x = referenceTable_split,
            y = referenceDictionary)

referenceTable_split <- referenceTable_split %>%
  select(referenceOriginal_split = referenceOriginal)

if (nrow(referenceTable_split) == 0)
  referenceTable_split <-
  rbind(referenceTable_split, list(NA))

### original
referenceTable_original <- dbTable %>%
  filter(is.na(referenceOriginal_doi)) %>%
  filter(is.na(referenceOriginal_pubmed)) %>%
  filter(is.na(referenceOriginal_title)) %>%
  filter(is.na(referenceOriginal_publishingDetails)) %>%
  filter(is.na(referenceOriginal_split)) %>%
  filter(!is.na(referenceOriginal_original)) %>%
  distinct(referenceOriginal_original) %>%
  select(referenceOriginal = referenceOriginal_original)

if (file.exists(pathDataInterimDictionariesReferenceDictionary))
  referenceTable_original <-
  anti_join(x = referenceTable_original,
            y = referenceDictionary)

referenceTable_original <- referenceTable_original %>%
  select(referenceOriginal_original = referenceOriginal)

if (nrow(referenceTable_original) == 0)
  referenceTable_original <-
  rbind(referenceTable_original, list(NA))

### full
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
  select(organismOriginal,
         referenceType,
         referenceValue)

## full table
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

# exporting
print(x = "exporting")
## creating directories if they do not exist
### interim tables
ifelse(
  test = !dir.exists(pathDataInterimTables),
  yes = dir.create(pathDataInterimTables),
  no = FALSE
)

#### original
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginal),
  yes = dir.create(pathDataInterimTablesOriginal),
  no = FALSE
)

##### organism
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalOrganism),
  yes = dir.create(pathDataInterimTablesOriginalOrganism),
  no = file.remove(
    list.files(path = pathDataInterimTablesOriginalOrganism,
               full.names = TRUE)
  ) &
    dir.create(pathDataInterimTablesOriginalOrganism,
               showWarnings = FALSE)
)

##### reference
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReference),
  yes = dir.create(pathDataInterimTablesOriginalReference),
  no = FALSE
)

##### reference folders
###### original
ifelse(
  test = !dir.exists(pathDataInterimTablesOriginalReferenceOriginalFolder),
  yes = dir.create(pathDataInterimTablesOriginalReferenceOriginalFolder),
  no = file.remove(
    list.files(path = pathDataInterimTablesOriginalReferenceOriginalFolder,
               full.names = TRUE)
  )  &
    dir.create(pathDataInterimTablesOriginalOrganism,
               showWarnings = FALSE)
)

###### title
ifelse(
  !dir.exists(pathDataInterimTablesOriginalReferenceTitleFolder),
  dir.create(pathDataInterimTablesOriginalReferenceTitleFolder),
  no = file.remove(
    list.files(path = pathDataInterimTablesOriginalReferenceTitleFolder,
               full.names = TRUE)
  )  &
    dir.create(pathDataInterimTablesOriginalOrganism,
               showWarnings = FALSE)
)

##### structure
ifelse(
  !dir.exists(pathDataInterimTablesOriginalStructure),
  dir.create(pathDataInterimTablesOriginalStructure),
  FALSE
)

## writing files
### organism
if (nrow(organismTable) != 0)
  split_data_table(
    x = organismTable,
    no_rows_per_frame = 10000,
    text = "",
    path_to_store = pathDataInterimTablesOriginalOrganism
  )

### reference
#### DOI
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

#### pubmed
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

#### title
split_data_table(
  x = referenceTable_title,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceTitleFolder
)

#### publishingDetails
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

#### split
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

#### original
split_data_table(
  x = referenceTable_original,
  no_rows_per_frame = 1000,
  text = "",
  path_to_store = pathDataInterimTablesOriginalReferenceOriginalFolder
)

#### full
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

### structure
#### inchi
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

#### nominal
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

#### smiles
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

#### full
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

## table
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
