source("r/log_debug.R")
log_debug("This script checks for organism presence in cleaned reference title")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(data.table)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

log_debug("... functions")
source("r/y_as_na.R")

log_debug("loading crossref translations file, this may take a while")
dataTranslated <-
  read_delim(
    file = pathDataInterimTablesTranslatedReferenceFile,
    delim = "\t",
    col_types = cols(.default = "c"),
    col_select =
      c(
        "organismType",
        "organismValue",
        "referenceType",
        "referenceValue",
        "organismDetected",
        "referenceCleanedType" = "referenceTranslatedType",
        "referenceCleanedValue" = "referenceTranslatedValue",
        "level"
      )
  ) %>%
  data.frame()
dataTranslated
log_debug("checking for organism in title, may take a while if running full mode")
dataCleanedScore <- dataTranslated %>%
  filter(referenceCleanedType == "title") %>%
  filter(!is.na(organismValue) &
    !is.na(referenceCleanedValue)) %>%
  mutate(referenceCleaned_scoreTitleOrganism = ifelse(
    test = str_detect(
      string = fixed(tolower(
        referenceCleanedValue
      )),
      pattern = fixed(tolower(
        word(organismDetected, 1)
      ))
    ),
    yes = 1,
    no = 0
  )) %>%
  pivot_wider(
    names_from = referenceCleanedType,
    names_prefix = "referenceCleaned_",
    values_from = referenceCleanedValue
  ) %>%
  mutate_all(as.character) %>%
  distinct() 
  dataCleanedScore
  dataCleanedScore <- dataCleanedScore %>% 
    pivot_longer(
    cols = starts_with("referenceCleaned_"),
    names_to = c("drop", "referenceCleanedType"),
    names_sep = "_",
    values_to = "referenceCleanedValue",
    values_drop_na = TRUE
  ) %>%
  filter(referenceCleanedType == "scoreTitleOrganism") %>%
  filter(referenceCleanedValue == 1) %>%
  select(-drop) %>%
  distinct()

dataCleanedJoined <- bind_rows(dataTranslated, dataCleanedScore) %>%
  filter(!is.na(organismValue))

rm(dataTranslated)

# this is because sadly crossref does not always give the same score, therefore
## we do not have unique values ...
subDataCleanedJoined_1 <- dataCleanedJoined %>%
  filter(referenceCleanedType == "scoreCrossref") %>%
  distinct(
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    organismDetected,
    level,
    .keep_all = TRUE
  )

# this is because sadly crossref does not always give the same DOI, therefore
## we do not have unique values ...
subDataCleanedJoined_2 <- dataCleanedJoined %>%
  filter(referenceCleanedType != "scoreCrossref") %>%
  distinct(
    organismType,
    organismValue,
    referenceType,
    referenceValue,
    organismDetected,
    level,
    referenceCleanedType,
    .keep_all = TRUE
  )

rm(dataCleanedJoined)

dataCleanedJoinedUnique <-
  bind_rows(
    subDataCleanedJoined_1,
    subDataCleanedJoined_2
  )

rm(
  subDataCleanedJoined_1,
  subDataCleanedJoined_2
)

gc()

log_debug("manipulating and keeping best result only (long step)")
dataCleanedJoinedWide <- dataCleanedJoinedUnique %>%
  pivot_wider(
    names_from = referenceCleanedType,
    names_prefix = "referenceCleaned_",
    values_from = referenceCleanedValue
  )

rm(dataCleanedJoinedUnique)

dataCleanedJoinedWide <- dataCleanedJoinedWide %>%
  union_all(data.frame(referenceCleaned_scoreTitleOrganism = character()))

dataCleanedJoinedWide_1 <- dataCleanedJoinedWide %>%
  filter(referenceType == "doi" |
    referenceType == "pubmed" |
    referenceType == "title") %>%
  group_by(organismType, organismValue, organismDetected, referenceValue) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreCrossref))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreTitleOrganism))) %>%
  arrange(as.numeric(referenceCleaned_scoreDistance)) %>%
  ungroup() %>%
  distinct(organismType,
    organismValue,
    organismDetected,
    referenceValue,
    .keep_all = TRUE
  ) %>%
  select(-level) %>%
  mutate(
    referenceCleaned_scoreComplement_date = NA,
    referenceCleaned_scoreComplement_author = NA,
    referenceCleaned_scoreComplement_journal = NA,
    referenceCleaned_scoreComplement_total = NA
  )

dataCleanedJoinedWide_2 <- dataCleanedJoinedWide %>%
  filter(
    referenceType == "original" |
      referenceType == "publishingDetails" |
      referenceType == "split"
  ) %>%
  mutate(
    referenceCleaned_scoreComplement_date = ifelse(
      str_detect(
        string = referenceValue,
        pattern = substr(x = referenceCleaned_date, start = 1, stop = 4)
      ),
      yes = 1,
      no = 0
    ),
    referenceCleaned_scoreComplement_author = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(word(
          referenceCleaned_author, 1
        )))
      ),
      yes = 1,
      no = 0
    ),
    referenceCleaned_scoreComplement_journal = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(referenceCleaned_journal))
      ),
      yes = 1,
      no = 0
    )
  ) %>%
  mutate(
    referenceCleaned_scoreComplement_total =
      referenceCleaned_scoreComplement_date +
        referenceCleaned_scoreComplement_author +
        referenceCleaned_scoreComplement_journal
  ) %>%
  group_by(
    organismType,
    organismValue,
    organismDetected,
    referenceValue
  ) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreCrossref))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreComplement_total))) %>%
  arrange(desc(as.numeric(referenceCleaned_scoreTitleOrganism))) %>%
  ungroup() %>%
  distinct(organismType,
    organismValue,
    organismDetected,
    referenceValue,
    .keep_all = TRUE
  ) %>%
  select(-level)

rm(dataCleanedJoinedWide)

dataCleanedJoinedWideScore <- bind_rows(
  dataCleanedJoinedWide_1,
  dataCleanedJoinedWide_2
) %>%
  mutate(referenceCleaned_doi = toupper(referenceCleaned_doi))

rm(
  dataCleanedJoinedWide_1,
  dataCleanedJoinedWide_2
)

subDataClean_doi <- dataCleanedJoinedWideScore %>%
  filter(!is.na(referenceCleaned_doi)) %>%
  distinct(referenceCleaned_doi) %>%
  mutate_all(as.character)

subDataClean_pmid <- dataCleanedJoinedWideScore %>%
  filter(referenceType == "pubmed") %>%
  distinct(referenceValue) %>%
  mutate_all(as.character)

if (mode == "test") {
  PMC_ids <-
    data.frame(
      DOI = NA,
      PMCID = NA,
      PMID = NA
    ) %>%
    mutate_all(as.character) %>%
    tibble()
} else if (mode == "manual") {
  PMC_ids <-
    data.frame(
      DOI = NA,
      PMCID = NA,
      PMID = NA
    ) %>%
    mutate_all(as.character) %>%
    tibble()
} else {
  log_debug("loading pmcid file, this may take a while")
  # here because of memory
  PMC_ids <- read_delim(
    file = pathDataExternalTranslationSourcePubmedFile,
    delim = ",",
    col_types = cols(.default = "c")
  ) %>%
    filter(!is.na(DOI) | !is.na(PMID)) %>%
    select(
      DOI,
      PMCID,
      PMID
    ) %>%
    mutate(DOI = toupper(DOI)) %>%
    mutate_all(as.character) %>%
    tibble()
}

log_debug("adding PMID and PMCID")
df_doi <- left_join(subDataClean_doi,
  PMC_ids,
  by = c("referenceCleaned_doi" = "DOI")
) %>%
  filter(!is.na(PMID) | !is.na(PMCID)) %>%
  select(
    referenceCleaned_doi,
    referenceCleaned_pmid = PMID,
    referenceCleaned_pmcid = PMCID
  )

df_pubmed <- left_join(subDataClean_pmid,
  PMC_ids,
  by = c("referenceValue" = "PMID")
) %>%
  filter(!is.na(referenceValue) | !is.na(PMCID)) %>%
  select(referenceValue,
    referenceCleaned_pmcid = PMCID
  ) %>%
  mutate(referenceCleaned_pmid = referenceValue)

tableJoined <- left_join(dataCleanedJoinedWideScore, df_doi)

referenceTable <-
  left_join(tableJoined,
    df_pubmed,
    by = c("referenceValue" = "referenceValue")
  ) %>%
  mutate(
    referenceCleaned_pmid = ifelse(
      test = !is.na(referenceCleaned_pmid.x),
      yes = referenceCleaned_pmid.x,
      no = referenceCleaned_pmid.y
    ),
    referenceCleaned_pmcid = ifelse(
      test = !is.na(referenceCleaned_pmcid.x),
      yes = referenceCleaned_pmcid.x,
      no = referenceCleaned_pmcid.y
    )
  ) %>%
  select(
    organismType,
    organismValue,
    organismDetected,
    referenceType,
    referenceValue,
    referenceCleanedDoi = referenceCleaned_doi,
    referenceCleanedPmcid = referenceCleaned_pmcid,
    referenceCleanedPmid = referenceCleaned_pmid,
    referenceCleanedTitle = referenceCleaned_title,
    referenceCleaned_journal,
    referenceCleaned_date,
    referenceCleaned_author,
    referenceCleaned_score_crossref = referenceCleaned_scoreCrossref,
    referenceCleaned_score_distance = referenceCleaned_scoreDistance,
    referenceCleaned_score_titleOrganism = referenceCleaned_scoreTitleOrganism,
    referenceCleaned_score_complementDate = referenceCleaned_scoreComplement_date,
    referenceCleaned_score_complementAuthor = referenceCleaned_scoreComplement_author,
    referenceCleaned_score_complementJournal = referenceCleaned_scoreComplement_journal,
    referenceCleaned_score_complementTotal = referenceCleaned_scoreComplement_total
  ) %>%
  distinct() %>%
  mutate(across(everything(), ~ y_as_na(.x, "NULL")))

log_debug("ensuring directories exist")
ifelse(
  test = !dir.exists(pathDataInterimTablesProcessed),
  yes = dir.create(pathDataInterimTablesProcessed),
  no = paste(pathDataInterimTablesProcessed, "exists")
)

ifelse(
  test = !dir.exists(pathDataInterimTablesProcessedReference),
  yes = dir.create(pathDataInterimTablesProcessedReference),
  no = paste(pathDataInterimTablesProcessedReference, "exists")
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedReferenceFile)
write_delim(
  x = referenceTable,
  file = pathDataInterimTablesProcessedReferenceFile,
  delim = "\t"
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
