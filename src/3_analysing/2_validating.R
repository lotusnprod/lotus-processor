## this is very dirty for now, I'll be cleaning it later on

cat(
  "This script aims to establish filtering criteria to validate \n",
  "documented pairs according to manually analyzed ones \n"
)

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(stringi)
library(tidyverse)
library(plotly)

cat("... functions \n")
source("r/filter.R")
source("r/myDirtyValidationFig.R")
source("r/vroom_safe.R")

cat("loading files ... \n")
oldDbNames <-
  read_delim(
    file = "../data/interim/dictionaries/dbNames.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_AR_old <-
  read_delim(
    file = "../data/validation/old/AR.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "AR")

sampleAllONPDB_PMA_old <-
  read_delim(
    file = "../data/validation/old/PMA.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "PMA")

sampleAllONPDB_AR <-
  read_delim(
    file = "../data/validation/new/AR.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "AR")

sampleAllONPDB_JB <-
  read_delim(
    file = "../data/validation/new/JB.csv",
    delim = ",",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "JB") %>%
  select(1:20)

sampleAllONPDB_PMA <-
  read_delim(
    file = "../data/validation/new/PM.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "PMA")

sampleAllONPDB_publishingDetails <-
  read_delim(
    file = "../data/validation/new/publishingDetails.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_additionalSet <-
  read_delim(
    file = "../data/validation/new/additionalSet.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_additionalSetBis <-
  read_delim(
    file = "../data/validation/new/additionalSetBis.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleCondifent_PMA <-
  read_delim(
    file = "../data/validation/confident/100confidentPMAChecked.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(curator == "PMA") %>%
  mutate(curator = "PMA2")

cat("... documented pairs \n")
inhouseDbMinimal <-
  vroom_read_safe(path = pathDataInterimTablesCuratedTable)

cat("... reference metadata \n")
structureMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesStructureMetadata) %>%
  distinct(
    structureCleanedSmiles,
    structureCleaned_smiles2D,
    structureCleanedInchi,
    structureCleaned_inchi2D,
    structureCleanedInchikey,
    structureCleaned_inchikey2D,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  )

cat("... organism metadata \n")
organismMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismMetadata) %>%
  distinct(
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

cat("... reference metadata \n")
referenceMetadata <-
  vroom_read_safe(path = pathDataInterimDictionariesReferenceMetadata) %>%
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
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism,
    referenceCleaned_score_complementTotal
  )

sampleAllONPDB_old <-
  rbind(sampleAllONPDB_AR_old, sampleAllONPDB_PMA_old)

colnames(sampleAllONPDB_old)[3] <- "structure_inchi"
colnames(sampleAllONPDB_old)[4] <- "structure_smiles"
colnames(sampleAllONPDB_old)[5] <- "structure_nominal"
colnames(sampleAllONPDB_old)[6] <- "reference_doi"
colnames(sampleAllONPDB_old)[7] <- "reference_external"
colnames(sampleAllONPDB_old)[8] <- "reference_pubmed"
colnames(sampleAllONPDB_old)[9] <- "reference_title"
colnames(sampleAllONPDB_old)[10] <- "reference_original"

table_old <- sampleAllONPDB_old %>%
  mutate_all(as.character) %>%
  mutate(
    validated = gsub(
      pattern = "NOT_OK",
      replacement = "N",
      x = validated
    ),
    validated = gsub(
      pattern = "OK",
      replacement = "Y",
      x = validated
    )
  ) %>%
  pivot_longer(
    6:10,
    names_to = c("origin1", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    3:5,
    names_to = c("origin2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) %>%
  filter(!is.na(referenceValue)) %>%
  group_by(inchikeySanitized) %>%
  distinct(referenceCleanedDoi, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-origin1, -origin2) %>%
  select(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned = organismLowestTaxon,
    structureCleanedInchi = inchiSanitized,
    structureCleanedInchikey = inchikeySanitized,
    structureCleanedSmiles = smilesSanitized,
    referenceCleanedDoi,
    referenceCleanedTitle,
    curator,
    validated,
    comments
  )

sampleAllONPDB <- bind_rows(
  sampleAllONPDB_AR,
  sampleAllONPDB_JB,
  sampleAllONPDB_PMA,
  sampleAllONPDB_publishingDetails,
  sampleAllONPDB_additionalSet,
  sampleAllONPDB_additionalSetBis,
  sampleCondifent_PMA
)

table <- sampleAllONPDB %>%
  mutate_all(as.character) %>%
  mutate(
    validated = gsub(
      pattern = "NotOK",
      replacement = "N",
      x = validated
    ),
    validated = gsub(
      pattern = "OK",
      replacement = "Y",
      x = validated
    ),
    validated = gsub(
      pattern = "NO",
      replacement = "N",
      x = validated
    ),
    validated = gsub(
      pattern = "YES",
      replacement = "Y",
      x = validated
    ),
    validated = gsub(
      pattern = "ToCheck",
      replacement = "Y/N",
      x = validated
    ),
    validated = gsub(
      pattern = "HUM",
      replacement = "Y/N",
      x = validated
    )
  ) %>%
  select(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey = structureCleanedInchikey3D,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    curator,
    validated,
    comments
  )

globalSample <- bind_rows(table_old, table) %>%
  mutate(organismValue = organismOriginal) %>%
  distinct(
    database,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleanedInchikey,
    referenceCleanedDoi,
    curator,
    validated,
    comments
  ) %>%
  filter(!is.na(validated)) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi)) %>%
  left_join(., inhouseDbMinimal) %>%
  select(
    -referenceCleanedTitle,
    -organismType
  ) %>%
  distinct()

a <- paste0("\\b", oldDbNames$oldDbName, "\\b")
b <- oldDbNames$newDbName

globalSample$database <- stri_replace_all_regex(
  str = globalSample$database,
  pattern = a,
  replacement = b,
  case_insensitive = FALSE,
  vectorize_all = FALSE
)

cat("adding metadata \n")
inhouseDbFull <- inhouseDbMinimal %>%
  left_join(., structureMetadata) %>%
  left_join(., organismMetadata) %>%
  left_join(., referenceMetadata)

cat("joining manual validation results with documented pairs \n")
realMetaSample <- inner_join(globalSample, inhouseDbFull) %>%
  distinct(
    database,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    curator,
    validated,
    comments,
    .keep_all = TRUE
  )

cat("filtering results ... \n")
cat("... validated set \n")
realSampleFiltered <-
  filter_dirty(dataframe = realMetaSample)

cat("... rejected set \n")
antiFilter <- anti_join(realMetaSample, realSampleFiltered)

cat("counting results ... \n")
cat("... per category on validation set \n")
table_count <- myDirtyF(table = realMetaSample)

cat("... global on validation set \n")
table_count_global <- myDirtyC(table = realMetaSample)

cat("... per category on validated set \n")
tableFiltered_count <-
  myDirtyF(table = realSampleFiltered)

cat("... global on validated set \n")
tableFiltered_count_global <-
  myDirtyC(table = realSampleFiltered)

cat("... per category on rejected set \n")
tableAntiFiltered_count <-
  myDirtyF(table = antiFilter)

cat("... global on rejected set \n")
tableAntiFiltered_count_global <-
  myDirtyC(table = antiFilter)

cat("visualizing ... \n")
cat("... validation set per category \n")
fig_full <-
  myDirtyP(
    table = table_count,
    yaxismax = 140,
    title = "full version"
  )
fig_full

cat("... validated set per category \n")
fig_filtered <-
  myDirtyP(
    table = tableFiltered_count,
    yaxismax = 140,
    title = "filtered version"
  )
fig_filtered

cat("... rejected set per category \n")
fig_anti <-
  myDirtyP(
    table = tableAntiFiltered_count,
    yaxismax = 140,
    title = "anti version"
  )
fig_anti

cat("... validation set global \n")
newfull <- myDirtyQ(
  table = table_count_global,
  yaxismax = 550,
  title = "new full version"
)
newfull

cat("... validated set global \n")
newfiltered <- myDirtyQ(
  table = tableFiltered_count_global,
  yaxismax = 550,
  title = "new filtered version"
)
newfiltered

cat("... rejected set global \n")
antifull <- myDirtyQ(
  table = tableAntiFiltered_count_global,
  yaxismax = 550,
  title = "anti full version"
)
antifull

cat("calculating statistics ... \n")
old <- table_count %>%
  select(referenceType,
    tot,
    y,
    n,
    ratio1 = ratio
  )

cat("... true positives and false positives \n")
new <- tableFiltered_count %>%
  select(
    referenceType,
    fil = tot,
    tp = y,
    fp = n,
    ratio2 = ratio
  )

cat("... true negatives and false negatives \n")
anti <- tableAntiFiltered_count %>%
  select(
    referenceType,
    anti = tot,
    fn = y,
    tn = n,
    ratioAnti = ratio
  )

f1Table <- full_join(old, new) %>%
  replace(is.na(.), 0)

beta <- 0.5

cat("... precision, recall and Fbeta", beta, "score \n")
f1Table <- full_join(f1Table, anti) %>%
  mutate(
    tpfn = tp + fn,
    tpfp = tp + fp
  ) %>%
  mutate(
    recall = tp / tpfn,
    precision = tp / tpfp
  ) %>%
  mutate(
    rxp = recall * precision,
    rpp = recall + precision
  ) %>%
  mutate(f2 = 2 * rxp / rpp) %>%
  select(-tpfn, -tpfp, -rxp, -rpp) %>%
  mutate(fbeta = (1 + beta^2) * (precision * recall) / ((precision * beta^
    2) + recall)) %>%
  replace(. == "NaN", 0)


cat("... correcting with references ratios \n")
refRatios <- inhouseDbFull %>%
  filter(database != "dnp") %>%
  group_by(referenceType) %>%
  count() %>%
  mutate(prop = n / sum(.$n)) %>%
  select(-n)

f2Table <- left_join(f1Table, refRatios) %>%
  mutate(
    f2corr_temp = f2 * prop,
    recall_temp = recall * prop,
    precision_temp = precision * prop
  ) %>%
  mutate(
    f2corr = sum(.$f2corr_temp),
    recallcorr = sum(.$recall_temp),
    precisioncorr = sum(.$precision_temp)
  ) %>%
  mutate()

f2Table_reworked <- f2Table %>%
  arrange(desc(prop)) %>%
  select(
    `reference type` = referenceType,
    `true positives` = tp,
    `false positives` = fp,
    `false negatives` = fn,
    `true negatives` = tn,
    `relative abundance` = prop,
    precision,
    recall,
    `F0.5 score` = f2,
    recallcorr,
    precisioncorr,
    f2corr
  ) %>%
  rbind(
    .,
    tibble(
      `reference type` = "Total",
      `true positives` = sum(.$`true positives`),
      `false positives` = sum(.$`false positives`),
      `false negatives` = sum(.$`false negatives`),
      `true negatives` = sum(.$`true negatives`),
      `relative abundance` = sum(.$`relative abundance`),
    ),
    tibble(
      `reference type` = "Corrected total",
      `precision` = unique(.$precisioncorr),
      `recall` = unique(.$recallcorr),
      `F0.5 score` = unique(.$f2corr),
    )
  ) %>%
  select(
    `reference type`,
    `true positives`,
    `false positives`,
    `false negatives`,
    `true negatives`,
    `relative abundance`,
    precision,
    recall,
    `F0.5 score`,
  )

cat("applying the filtering criteria to the whole DB, this may take a while \n")
openDb <- inhouseDbFull %>%
  filter_dirty() %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleaned_inchi2D,
    structureCleaned_smiles2D,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  mutate(structureCleaned_inchikey2D = substring(
    text = structureCleanedInchikey,
    first = 1,
    last = 14
  )) %>%
  select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleaned_inchikey2D,
    structureCleaned_inchi2D,
    structureCleaned_smiles2D,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

cat("outputting dnp pairs")
dnpDb <- inhouseDbFull %>%
  filter(database == "dnp") %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleaned_inchi2D,
    structureCleaned_smiles2D,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  mutate(structureCleaned_inchikey2D = substring(
    text = structureCleanedInchikey,
    first = 1,
    last = 14
  )) %>%
  select(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleaned_inchikey2D,
    structureCleaned_inchi2D,
    structureCleaned_smiles2D,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

cat("outputing correct entries from manually validated set \n")
manuallyValidatedSet <- realMetaSample %>%
  filter(validated == "Y")

cat("outputing incorrect entries from validated set \n")
manuallyRemovedEntries <- realMetaSample %>%
  filter(validated != "Y")

openDbClean <- anti_join(openDb, manuallyRemovedEntries)

if (mode != "test") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  validationSet <- anti_join(openDbClean, realMetaSample) %>%
    sample_n(100)
}

cat("loading validation set \n")
validationSetFilled_1 <-
  read_delim(
    file = "../data/validation/validationSet.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(!is.na(validated)) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

if (mode != "test") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  validationSet2 <-
    anti_join(openDbClean, validationSetFilled_1) %>%
    sample_n(13)
}

cat("loading validation set bis \n")
validationSetFilled_2 <-
  read_delim(
    file = "../data/validation/validationSetBis.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(!is.na(validated)) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

if (mode != "test") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  validationSet3 <-
    anti_join(openDbClean, validationSetFilled_1) %>%
    anti_join(., validationSetFilled_2) %>%
    sample_n(19)
}

cat("loading validation set ter \n")
validationSetFilled_3 <-
  read_delim(
    file = "../data/validation/validationSetTer.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  filter(!is.na(validated)) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

validationSetFilled <-
  bind_rows(
    validationSetFilled_1,
    validationSetFilled_2,
    validationSetFilled_3
  ) %>%
  mutate(organismValue = organismOriginal) %>%
  select(
    structureCleanedInchikey = structureCleanedInchikey3D,
    everything()
  ) %>%
  distinct(
    database,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    curator,
    validated,
    comments
  ) %>%
  left_join(., inhouseDbMinimal) %>%
  select(-referenceCleanedTitle, -organismType) %>%
  distinct()

validationSetFilled$database <- stri_replace_all_regex(
  str = validationSetFilled$database,
  pattern = a,
  replacement = b,
  case_insensitive = FALSE,
  vectorize_all = FALSE
)

realValidationSetFilled <-
  inner_join(validationSetFilled, openDbClean) %>%
  distinct(
    database,
    organismType,
    organismValue,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    curator,
    validated,
    comments,
    .keep_all = TRUE
  )

finalStats <- realValidationSetFilled %>%
  group_by(referenceType) %>%
  count(validated == "Y") %>%
  ungroup()

finalStats_reworked <- finalStats %>%
  mutate(
    `true positives` = ifelse(
      test = `validated == "Y"` == TRUE,
      yes = n,
      no = 0
    ),
    `false positives` = ifelse(
      test = `validated == "Y"` == FALSE,
      yes = n,
      no = 0
    )
  ) %>%
  select(
    `reference type` = referenceType,
    everything(),
    -n,
    -`validated == "Y"`
  ) %>%
  group_by(`reference type`) %>%
  summarise(
    `true positives Validation` = max(`true positives`),
    `false positives Validation` = max(`false positives`)
  ) %>%
  rbind(
    .,
    tibble(
      `reference type` = "Total",
      `true positives Validation` = sum(.$`true positives Validation`),
      `false positives Validation` = sum(.$`false positives Validation`)
    )
  )

finalTable <- left_join(f2Table_reworked, finalStats_reworked) %>%
  mutate_if(.predicate = is.numeric, ~ round(., digits = 2))

cat("outputing correct entries from manually validated set \n")
manuallyValidatedSet2 <- realValidationSetFilled %>%
  filter(validated == "Y")

manuallyValidatedSet3 <-
  bind_rows(manuallyValidatedSet, manuallyValidatedSet2)

cat("outputing incorrect entries from validated set \n")
manuallyRemovedEntries2 <- realValidationSetFilled %>%
  filter(validated != "Y")

manuallyRemovedEntries3 <-
  bind_rows(manuallyRemovedEntries, manuallyRemovedEntries2)

openDbClean2 <- anti_join(openDbClean, manuallyRemovedEntries3) %>%
  filter(!database %in% forbidden_export)

cat("exporting \n")
if (mode == "full") {
  cat("../data/validation/manuallyValidated.tsv.gz", "\n")
}
if (mode == "full") {
  write.table(
    x = manuallyValidatedSet3,
    file = gzfile(
      description = "../data/validation/manuallyValidated.tsv.gz",
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (mode == "full") {
  cat("../data/validation/manuallyRemoved.tsv.gz", "\n")
}
if (mode == "full") {
  write.table(
    x = manuallyRemovedEntries3,
    file = gzfile(
      description = "../data/validation/manuallyRemoved.tsv.gz",
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )

  write.table(
    x = finalTable,
    file = "../data/validation/tableStats.csv",
    row.names = FALSE,
    quote = FALSE,
    sep = ",",
    fileEncoding = "UTF-8"
  )
}

cat(pathDataInterimTablesAnalysedPlatinum, "\n")
vroom_write_safe(
  x = openDbClean2,
  path = pathDataInterimTablesAnalysedPlatinum
)

cat(file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz"), "\n")
vroom_write_safe(
  x = dnpDb,
  path = file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz")
)

cat(
  file.path(
    pathDataInterimTablesAnalysed,
    "validationSet.tsv"
  ),
  "\n"
)
if (exists("validationSet")) {
  write.table(
    x = validationSet,
    file = file.path(
      pathDataInterimTablesAnalysed,
      "validationSet.tsv"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

cat(
  file.path(
    pathDataInterimTablesAnalysed,
    "validationSetBis.tsv"
  ),
  "\n"
)

if (exists("validationSet2")) {
  write.table(
    x = validationSet2,
    file = file.path(
      pathDataInterimTablesAnalysed,
      "validationSetBis.tsv"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("validationSet3")) {
  write.table(
    x = validationSet3,
    file = file.path(
      pathDataInterimTablesAnalysed,
      "validationSetTer.tsv"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
