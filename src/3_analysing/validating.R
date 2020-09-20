# this is very dirty for now, I'll be cleaning it later on

source("paths.R")
library(tidyverse)
library(plotly)

sampleAllONPDB_AR_old <-
  read_delim(file = "../data/validation/old/AR.tsv",
             delim = "\t") %>%
  filter(curator == "AR")

sampleAllONPDB_PMA_old <-
  read_delim(file = "../data/validation/old/PMA.tsv",
             delim = "\t") %>%
  filter(curator == "PMA")

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

sampleAllONPDB_AR <-
  read_delim(file = "../data/validation/new/AR.tsv",
             delim = "\t") %>%
  filter(curator == "AR")

sampleAllONPDB_JB <-
  read_delim(file = "../data/validation/new/JB.csv",
             delim = ";") %>%
  filter(curator == "JB") %>%
  select(1:20)

sampleAllONPDB_PMA <-
  read_delim(file = "../data/validation/new/PM.tsv",
             delim = "\t") %>%
  filter(curator == "PMA")

sampleAllONPDB_publishingDetails <-
  read_delim(file = "../data/validation/new/publishingDetails.tsv",
             delim = "\t")

sampleCondifent_PMA <-
  read_delim(file = "../data/validation/confident/100confidentPMAChecked.tsv",
             delim = "\t") %>%
  filter(curator == "PMA") %>%
  mutate(curator = "PMA2")

sampleAllONPDB <- bind_rows(
  sampleAllONPDB_AR,
  sampleAllONPDB_JB,
  sampleAllONPDB_PMA,
  sampleAllONPDB_publishingDetails,
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
         -referenceCleanedTitle)

openDb <- read_delim(
  file = gzfile(
    file.path(pathDataInterimTablesAnalysed, "openDbTriplets.tsv.gz")
  ),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
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
    referenceCleanedTitle
  ) %>%
  data.frame()

referenceMetadata <-
  read_delim(file = gzfile(pathDataInterimDictionariesReferenceMetadata),
             delim = "\t") %>%
  select(
    organismCleaned,
    referenceCleanedDoi,
    referenceCleaned_score_crossref,
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism,
    referenceCleaned_score_complementTotal,
    referenceCleaned_score_complementDate,
    referenceCleaned_score_complementAuthor,
    referenceCleaned_score_complementJournal
  ) %>%
  distinct(organismCleaned, referenceCleanedDoi, .keep_all = TRUE)

realSample <- inner_join(globalSample, openDb) %>%
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

realMetaSample <- left_join(realSample, referenceMetadata)

realSampleFilteredBioTitle <- realMetaSample %>%
  filter(
    referenceCleaned_score_crossref == 1 |
      referenceCleaned_score_distance <= 5 |
      referenceCleaned_score_complementTotal >= 2
  ) %>%
  filter(referenceCleaned_score_titleOrganism == 1)

antiFilter <- anti_join(realMetaSample, realSampleFilteredBioTitle)

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

table_count <- myDirtyF(table = realSample)

table_count_global <- myDirtyC(table = realSample)

tableFiltered_count <-
  myDirtyF(table = realSampleFilteredBioTitle)

tableFiltered_count_global <-
  myDirtyC(table = realSampleFilteredBioTitle %>% filter(curator == "PMA2"))

tableAntiFiltered_count <-
  myDirtyF(table = antiFilter)

tableAntiFiltered_count_global <-
  myDirtyC(table = antiFilter)

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

fig_full <-
  myDirtyP(table = table_count,
           yaxismax = 140,
           title = "full version")
fig_full

fig_filtered <-
  myDirtyP(table = tableFiltered_count,
           yaxismax = 120,
           title = "filtered version")
fig_filtered

fig_anti <-
  myDirtyP(table = tableAntiFiltered_count,
           yaxismax = 120,
           title = "anti version")
fig_anti

newfull <- myDirtyQ(table = table_count_global,
                    yaxismax = 350,
                    title = "new full version")
newfull

newfiltered <- myDirtyQ(table = tableFiltered_count_global,
                        yaxismax = 120,
                        title = "new filtered version")
newfiltered

antifull <- myDirtyQ(table = tableAntiFiltered_count_global,
                     yaxismax = 350,
                     title = "new full version")
antifull

old <- table_count %>%
  select(
    referenceType,
    n = tot,
    tp1 = y,
    fp1 = n,
    ambiguous1 = mix,
    ratio1 = ratio
  )

new <- tableFiltered_count %>%
  select(
    referenceType,
    fil = tot,
    tp2 = y,
    fp2 = n,
    ambiguous2 = mix,
    ratio2 = ratio
  )

anti <- tableAntiFiltered_count %>%
  select(
    referenceType,
    anti = tot,
    fn1 = y,
    tn1 = n,
    ambiguousAnti = mix,
    ratioAnti = ratio
  )

f1Table <- full_join(old, new)

f1Table <- full_join(f1Table, anti) %>%
  mutate(tpfn = tp2 + fn1,
         tpfp = tp2 + fp2) %>%
  mutate(recall = tp2 / tpfn,
         precision = tp2 / tpfp) %>%
  mutate(rxp = recall * precision,
         rpp = recall + precision) %>%
  mutate(f1 = 2 * rxp / rpp)
