# this is very dirty for now, I'll be cleaning it later on

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

sampleCondifent_PMA <-
  read_delim(file = "../data/validation/confident/100confidentPMAChecked.tsv",
             delim = "\t") %>%
  filter(curator == "PMA") %>%
  mutate(curator = "PMA2")

sampleAllONPDB <- bind_rows(sampleAllONPDB_AR,
                            sampleAllONPDB_JB,
                            sampleAllONPDB_PMA,
                            sampleCondifent_PMA)

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
  relocate(c(referenceCleanedPmcid, referenceCleanedPmid), .after = referenceCleanedDoi)

openDb <- read_delim(
  file = gzfile("../data/interim/tables/4_analysed/openDbTriplets.tsv.gz"),
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
  ) %>%
  data.frame()

referenceMetadata <-
  read_delim(
    file = gzfile("../data/interim/dictionaries/reference/metadata.tsv.gz"),
    delim = "\t"
  ) %>%
  select(
    organismCleaned,
    referenceCleanedDoi,
    referenceCleaned_score_crossref,
    referenceCleaned_score_distance,
    referenceCleaned_score_titleOrganism
  ) %>%
  distinct(organismCleaned, referenceCleanedDoi, .keep_all = TRUE)

realSample <- inner_join(globalSample, openDb)

realMetaSample <- left_join(realSample, referenceMetadata)

realSampleFilteredBioTitle <- realMetaSample %>%
  filter(referenceCleaned_score_distance <= 5 |
           referenceType != "title") %>%
  filter(referenceCleaned_score_titleOrganism == 1)

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

tableFiltered_count <- myDirtyF(table = realSampleFilteredBioTitle)

tableFiltered_count_global <-
  myDirtyC(table = realSampleFilteredBioTitle %>% filter(curator == "PMA2"))

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

newfull <- myDirtyQ(table = table_count_global,
                    yaxismax = 350,
                    title = "new full version")
newfull

newfiltered <- myDirtyQ(table = tableFiltered_count_global,
                        yaxismax = 120,
                        title = "new filtered version")
newfiltered
