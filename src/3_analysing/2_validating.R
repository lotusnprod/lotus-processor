# this is very dirty for now, I'll be cleaning it later on

cat(
  "This script aims to establish filtering criteria to validate \n",
  "documented pairs according to manually analyzed ones \n"
)

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions.R")
source("functions/filter.R")

myDirtyF <- function(table) {
  table_tot <- table %>%
    group_by(referenceType) %>%
    count(name = "tot")
  
  table_y <- table %>%
    filter(validated == "Y") %>%
    group_by(referenceType) %>%
    count(name = "y")
  
  table_n <- table %>%
    filter(validated == "N") %>%
    group_by(referenceType) %>%
    count(name =  "n")
  
  table_mix <- table %>%
    filter(validated == "Y/N") %>%
    group_by(referenceType) %>%
    count(name = "mix")
  
  table_full <- left_join(table_tot, table_y)
  table_full <- left_join(table_full, table_n)
  table_full <- left_join(table_full, table_mix) %>%
    replace(is.na(.), 0) %>%
    mutate(ratio =  y / tot)
  
  return (table_full)
}
myDirtyC <- function(table) {
  table_tot <- table %>%
    count(name = "tot")
  
  table_y <- table %>%
    filter(validated == "Y") %>%
    count(name = "y")
  
  table_n <- table %>%
    filter(validated == "N") %>%
    count(name =  "n")
  
  table_mix <- table %>%
    filter(validated == "Y/N") %>%
    count(name = "mix")
  
  table_full <-
    bind_cols(table_tot, table_y, table_n, table_mix)  %>%
    replace(is.na(.), 0) %>%
    mutate(ratio =  y / tot)
  
  return (table_full)
}
myDirtyP <- function(table, title, yaxismax) {
  fig <-
    plot_ly(
      data = table,
      x = ~ referenceType,
      y = ~ y,
      type = 'bar',
      name = 'correct',
      color = I('green')
    ) %>%
    add_trace(y = ~ mix,
              name = 'ambiguous',
              color = I("orange")) %>%
    add_trace(
      y = ~ n,
      name = 'uncorrect',
      color = I('red'),
      text = ~ round(x = ratio, digits = 2),
      textposition = 'outside',
      textfont = list(color = I('black'), size = 20)
    ) %>%
    layout(
      title = title,
      yaxis = list(title = 'Count',
                   range = c(0, yaxismax)),
      barmode = 'stack'
    )
  return(fig)
}
myDirtyQ <- function(table, title, yaxismax) {
  fig <-
    plot_ly(
      data = table,
      y = ~ y,
      type = 'bar',
      name = 'correct',
      color = I('green')
    ) %>%
    add_trace(y = ~ mix,
              name = 'ambiguous',
              color = I("orange")) %>%
    add_trace(
      y = ~ n,
      name = 'uncorrect',
      color = I('red'),
      text = ~ round(x = ratio, digits = 2),
      textposition = 'outside',
      textfont = list(color = I('black'), size = 20)
    ) %>%
    layout(
      title = title,
      yaxis = list(title = 'Count',
                   range = c(0, yaxismax)),
      barmode = 'stack'
    )
  return(fig)
}

cat("... libraries \n")
library(plotly)

cat("loading files ... \n")
sampleAllONPDB_AR_old <-
  read_delim(file = "../data/validation/old/AR.tsv",
             delim = "\t") %>%
  filter(curator == "AR")

sampleAllONPDB_PMA_old <-
  read_delim(file = "../data/validation/old/PMA.tsv",
             delim = "\t") %>%
  filter(curator == "PMA")

sampleAllONPDB_AR <-
  read_delim(file = "../data/validation/new/AR.tsv",
             delim = "\t") %>%
  filter(curator == "AR")

sampleAllONPDB_JB <-
  read_delim(file = "../data/validation/new/JB.csv",
             delim = ",") %>%
  filter(curator == "JB") %>%
  select(1:20)

sampleAllONPDB_PMA <-
  read_delim(file = "../data/validation/new/PM.tsv",
             delim = "\t") %>%
  filter(curator == "PMA")

sampleAllONPDB_publishingDetails <-
  read_delim(file = "../data/validation/new/publishingDetails.tsv",
             delim = "\t")

sampleAllONPDB_additionalSet <-
  read_delim(file = "../data/validation/new/additionalSet.tsv",
             delim = "\t")

sampleCondifent_PMA <-
  read_delim(file = "../data/validation/confident/100confidentPMAChecked.tsv",
             delim = "\t") %>%
  filter(curator == "PMA") %>%
  mutate(curator = "PMA2")

cat("... documented pairs \n")
inhouseDbMinimal <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

cat("... organism metadata \n")
organismMetadata <-
  read_delim(
    file = gzfile(pathDataInterimDictionariesOrganismMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  select(1:5)

cat("... reference metadata \n")
referenceMetadata <-
  read_delim(
    file = gzfile(pathDataInterimDictionariesReferenceMetadata),
    delim = "\t",
    col_types = cols(.default = "c"),
    escape_double = FALSE,
    trim_ws = TRUE
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
  filter(
    !is.na(referenceCleanedDoi) |
      !is.na(referenceCleanedPmcid) |
      !is.na(referenceCleanedPmid)
  ) %>%
  relocate(c(referenceCleanedPmcid, referenceCleanedPmid), .after = referenceCleanedDoi) %>%
  select(-referenceCleanedPmcid,
         -referenceCleanedPmid,
         -referenceCleanedTitle) %>%
  filter(!is.na(validated))

cat("adding metadata \n")
cat("... organisms \n")
inhouseDbFull <- left_join(inhouseDbMinimal, organismMetadata)
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
realSampleFilteredBioTitle <-
  filter_dirty(dataframe = realMetaSample)

cat("... rejected set \n")
antiFilter <- anti_join(realMetaSample, realSampleFilteredBioTitle)

cat("counting results ... \n")
cat("... per category on validation set \n")
table_count <- myDirtyF(table = realMetaSample)

cat("... global on validation set \n")
table_count_global <- myDirtyC(table = realMetaSample)

cat("... per category on validated set \n")
tableFiltered_count <-
  myDirtyF(table = realSampleFilteredBioTitle)

cat("... global on validated set \n")
tableFiltered_count_global <-
  myDirtyC(table = realSampleFilteredBioTitle)

cat("... per category on rejected set \n")
tableAntiFiltered_count <-
  myDirtyF(table = antiFilter)

cat("... global on rejected set \n")
tableAntiFiltered_count_global <-
  myDirtyC(table = antiFilter)

cat("visualizing ... \n")
cat("... validation set per category \n")
fig_full <-
  myDirtyP(table = table_count,
           yaxismax = 140,
           title = "full version")
fig_full

cat("... validated set per category \n")
fig_filtered <-
  myDirtyP(table = tableFiltered_count,
           yaxismax = 140,
           title = "filtered version")
fig_filtered

cat("... rejected set per category \n")
fig_anti <-
  myDirtyP(table = tableAntiFiltered_count,
           yaxismax = 140,
           title = "anti version")
fig_anti

cat("... validation set global \n")
newfull <- myDirtyQ(table = table_count_global,
                    yaxismax = 500,
                    title = "new full version")
newfull

cat("... validated set global \n")
newfiltered <- myDirtyQ(table = tableFiltered_count_global,
                        yaxismax = 500,
                        title = "new filtered version")
newfiltered

cat("... rejected set global \n")
antifull <- myDirtyQ(table = tableAntiFiltered_count_global,
                     yaxismax = 500,
                     title = "anti full version")
antifull

cat("calculating statistics ... \n")
old <- table_count %>%
  select(referenceType,
         tot,
         y,
         n,
         ratio1 = ratio)

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
  mutate(tpfn = tp + fn,
         tpfp = tp + fp) %>%
  mutate(recall = tp / tpfn,
         precision = tp / tpfp) %>%
  mutate(rxp = recall * precision,
         rpp = recall + precision) %>%
  mutate(f2 = 2 * rxp / rpp) %>%
  select(-tpfn, -tpfp, -rxp, -rpp) %>%
  mutate(fbeta = (1 + beta ^ 2) * (precision * recall) / ((precision * beta ^
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
  mutate(structureCleanedInchikey2D = substring(text = structureCleanedInchikey3D,
                                                first = 1,
                                                last = 14)) %>%
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
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) %>%
  mutate(structureCleanedInchikey2D = substring(text = structureCleanedInchikey3D,
                                                first = 1,
                                                last = 14)) %>%
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
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

cat("outputing correct entries from manually validated set \n")
manuallyValidatedSet <- realSampleFilteredBioTitle %>%
  filter(validated == "Y")

cat("outputing incorrect entries from validated set \n")
manuallyRemovedEntries <- realSampleFilteredBioTitle %>%
  filter(validated != "Y")

openDbClean <- anti_join(openDb, manuallyRemovedEntries)

cat("exporting \n")
cat("../data/validation/manuallyValidated.tsv.gz", "\n")
write.table(
  x = manuallyValidatedSet,
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

cat("../data/validation/manuallyRemoved.tsv.gz", "\n")
write.table(
  x = manuallyRemovedEntries,
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

cat(pathDataInterimTablesAnalysedPlatinum, "\n")
write.table(
  x = openDbClean,
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

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
