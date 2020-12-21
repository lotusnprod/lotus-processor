# this is very dirty for now, I'll be cleaning it later on

cat(
  "This script aims to establish filtering criteria to validate \n",
  "documented pairs according to manually analyzed ones \n"
)

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(tidyverse)
library(plotly)

cat("... functions \n")
source("r/filter.R")
source("r/myDirtyValidationFig.R")

cat("loading files ... \n")
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
inhouseDbMinimal <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

cat("... reference metadata \n")
structureMetadata <-
  read_delim(
    file = gzfile(pathDataInterimDictionariesStructureMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  distinct(
    structureCleanedSmiles,
    structureCleanedInchi,
    structureCleanedInchikey3D,
    structureCleaned_inchikey2D,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  )

cat("... organism metadata \n")
organismMetadata <-
  read_delim(
    file = gzfile(pathDataInterimDictionariesOrganismMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  distinct(
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

cat("... reference metadata \n")
referenceMetadata <-
  read_delim(
    file = gzfile(pathDataInterimDictionariesReferenceMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
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
    structureCleanedInchikey3D = inchikeySanitized,
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
    structureCleanedInchikey3D,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    curator,
    validated,
    comments
  )

globalSample <- bind_rows(table_old, table) %>%
  distinct(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleanedInchikey3D,
    referenceCleanedDoi,
    curator,
    validated,
    comments
  ) %>%
  filter(!is.na(validated)) %>%
  mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

cat("adding metadata \n")
cat("... structures \n")
inhouseDbFull <-
  left_join(inhouseDbMinimal, structureMetadata)

cat("... organisms \n")
inhouseDbFull <-
  left_join(inhouseDbFull, organismMetadata)

cat("... references \n")
inhouseDbFull <- left_join(inhouseDbFull, referenceMetadata)

cat("joining manual validation results with documented pairs \n")
realMetaSample <- inner_join(globalSample, inhouseDbFull) %>%
  distinct(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey3D,
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
# β is chosen such that recall is considered β times as important as precision

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
  filter(database != "dnp_1") %>%
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

cat("applying the filtering criteria to the whole DB, this may take a while \n")
openDb <- inhouseDbFull %>%
  filter_dirty() %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  mutate(structureCleanedInchikey2D = substring(
    text = structureCleanedInchikey3D,
    first = 1,
    last = 14
  )) %>%
  select(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchikey2D,
    structureCleanedInchi,
    structureCleanedSmiles,
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
  filter(database == "dnp_1") %>%
  distinct(
    database,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchi,
    structureCleanedSmiles,
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
  mutate(structureCleanedInchikey2D = substring(
    text = structureCleanedInchikey3D,
    first = 1,
    last = 14
  )) %>%
  select(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    organismCleaned_dbTaxo,
    organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy,
    structureCleanedInchikey3D,
    structureCleanedInchikey2D,
    structureCleanedInchi,
    structureCleanedSmiles,
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

validationSetFilled <-
  bind_rows(validationSetFilled_1, validationSetFilled_2)

realValidationSetFilled <-
  inner_join(validationSetFilled, openDbClean) %>%
  distinct(
    database,
    organismOriginal,
    structureType,
    structureValue,
    referenceType,
    referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey3D,
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
  count(validated == "Y")

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
}

cat(pathDataInterimTablesAnalysedPlatinum, "\n")
write.table(
  x = openDbClean2,
  file = gzfile(
    description = pathDataInterimTablesAnalysedPlatinum,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz"), "\n")
write.table(
  x = dnpDb,
  file = gzfile(
    description = file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz"),
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
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

cat(file.path(
  pathDataInterimTablesAnalysed,
  "validationSetBis.tsv"
))
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
