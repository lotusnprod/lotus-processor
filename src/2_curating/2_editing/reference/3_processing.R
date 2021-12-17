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
dataCrossref <-
  read_delim(
    file = pathDataInterimDictionariesReferenceDictionary,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

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
  PMC_ids <- read_delim(
    file = pathDataExternalTranslationSourcePubmedFile,
    delim = ",",
    col_types = cols(.default = "c"),
    locale = locales
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

log_debug("loading references lists")
length <-
  length(list.files(
    path = pathDataInterimTablesTranslatedReference,
    pattern = "tsv$"
  ))

cut <- 100000

num <- as.integer(seq(
  from = 1 * cut,
  to = length * cut,
  by = cut
))

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

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesTranslatedReference,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      ".tsv"
    )

  outpath <-
    paste0(
      pathDataInterimTablesTranslatedReference,
      "/",
      str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      "_processed.tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  # log_debug("loading translated file")
  dataTranslated <-
    read_delim(
      file = inpath,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )

  dataCrossref_step <- dataCrossref %>%
    semi_join(
      dataTranslated,
      by = c(
        "referenceOriginal" = "referenceValue",
        "origin" = "referenceType"
      )
    ) %>%
    distinct(referenceOriginal,
      referenceTranslatedType,
      origin,
      level,
      .keep_all = TRUE
    ) %>%
    semi_join(
      dataTranslated,
      by = c(
        "referenceOriginal" = "referenceValue",
        "origin" = "referenceType"
      )
    ) %>%
    distinct(referenceOriginal,
      referenceTranslatedType,
      origin,
      level,
      .keep_all = TRUE
    )
  # log_debug("selecting minimal crossref data for DOI, PUBMED, TITLE")
  dataCrossref_1 <- dataCrossref_step %>%
    filter(
      referenceTranslatedType == "doi" |
        referenceTranslatedType == "title" |
        referenceTranslatedType == "scoreCrossref" |
        referenceTranslatedType == "scoreDistance"
    ) %>%
    pivot_wider(
      names_from = referenceTranslatedType,
      names_prefix = "referenceTranslated_",
      values_from = referenceTranslatedValue
    )

  # log_debug("selecting minimal crossref data for ORIGINAL, PUBDETAILS, SPLIT")
  dataCrossref_2 <- dataCrossref_step %>%
    filter(
      referenceTranslatedType == "doi" |
        referenceTranslatedType == "title" |
        referenceTranslatedType == "scoreCrossref" |
        referenceTranslatedType == "journal" |
        referenceTranslatedType == "date" |
        referenceTranslatedType == "author"
    ) %>%
    pivot_wider(
      names_from = referenceTranslatedType,
      names_prefix = "referenceTranslated_",
      values_from = referenceTranslatedValue
    )

  # log_debug("selecting minimal data for organism in title recognition")
  dataTranslated_title <- dataTranslated %>%
    filter(!is.na(organismDetected)) %>%
    inner_join(
      dataCrossref_step %>%
        filter(referenceTranslatedType == "title"),
      by = c(
        "referenceValue" = "referenceOriginal",
        "referenceType" = "origin"
      )
    )

  dataTranslated_title_distinct <- dataTranslated_title %>%
    distinct(organismDetected, referenceTranslatedValue)

  # log_debug("checking for organism in title, may take a while if running full mode")
  dataTranslated_title_distinct <- dataTranslated_title_distinct %>%
    mutate(referenceTranslated_scoreTitleOrganism = ifelse(
      test = str_detect(
        string = fixed(tolower(
          referenceTranslatedValue
        )),
        pattern = fixed(tolower(organismDetected))
      ),
      yes = 1,
      no = 0
    )) %>%
    filter(referenceTranslated_scoreTitleOrganism == 1) %>%
    distinct(
      organismDetected,
      referenceTranslated_title = referenceTranslatedValue,
      referenceTranslated_scoreTitleOrganism
    )

  # log_debug("splitting cases")
  dataTranslated_1 <- dataTranslated %>%
    filter(referenceType == "doi" |
      referenceType == "pubmed" |
      referenceType == "title") %>%
    left_join(
      dataCrossref_1,
      by = c(
        "referenceValue" = "referenceOriginal",
        "referenceType" = "origin"
      )
    )

  dataTranslated_2 <- dataTranslated %>%
    filter(
      referenceType == "original" |
        referenceType == "publishingDetails" |
        referenceType == "split"
    ) %>%
    left_join(
      dataCrossref_2,
      by = c(
        "referenceValue" = "referenceOriginal",
        "referenceType" = "origin"
      )
    )

  dataTranslated_3 <- dataTranslated %>%
    filter(
      referenceType != "doi" &
        referenceType != "pubmed" &
        referenceType != "title" &
        referenceType != "original" &
        referenceType != "publishingDetails" &
        referenceType != "split"
    )

  dataCleanedJoinedWide_1 <- dataTranslated_1 %>%
    distinct(
      organismType,
      organismValue,
      organismDetected,
      referenceType,
      referenceValue,
      referenceTranslated_title,
      referenceTranslated_scoreCrossref,
      referenceTranslated_scoreDistance,
      level
    ) %>%
    left_join(dataTranslated_title_distinct) %>%
    select(-referenceTranslated_title) %>%
    distinct() %>%
    group_by(
      organismType,
      organismValue,
      organismDetected,
      referenceType,
      referenceValue
    ) %>%
    arrange(desc(as.numeric(referenceTranslated_scoreCrossref))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreTitleOrganism
    ))) %>%
    arrange(as.numeric(referenceTranslated_scoreDistance)) %>%
    ungroup() %>%
    distinct(organismType,
      organismValue,
      organismDetected,
      referenceValue,
      .keep_all = TRUE
    ) %>%
    select(-level)

  dataCleanedJoinedWide_2_author <- dataTranslated_2 %>%
    distinct(
      referenceType,
      referenceValue,
      referenceTranslated_doi,
      referenceTranslated_title,
      referenceTranslated_author
    ) %>%
    mutate(referenceTranslated_scoreComplement_author = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(word(
          referenceTranslated_author, 1
        )))
      ),
      yes = 1,
      no = 0
    )) %>%
    filter(referenceTranslated_scoreComplement_author == 1)

  dataCleanedJoinedWide_2_date <- dataTranslated_2 %>%
    distinct(
      referenceType,
      referenceValue,
      referenceTranslated_doi,
      referenceTranslated_title,
      referenceTranslated_date
    ) %>%
    mutate(referenceTranslated_scoreComplement_date = ifelse(
      str_detect(
        string = referenceValue,
        pattern = substr(x = referenceTranslated_date, start = 1, stop = 4)
      ),
      yes = 1,
      no = 0
    )) %>%
    filter(referenceTranslated_scoreComplement_date == 1)

  dataCleanedJoinedWide_2_journal <- dataTranslated_2 %>%
    distinct(
      referenceType,
      referenceValue,
      referenceTranslated_doi,
      referenceTranslated_title,
      referenceTranslated_journal
    ) %>%
    mutate(referenceTranslated_scoreComplement_journal = ifelse(
      str_detect(
        string = tolower(referenceValue),
        pattern = fixed(tolower(referenceTranslated_journal))
      ),
      yes = 1,
      no = 0
    )) %>%
    filter(referenceTranslated_scoreComplement_journal == 1)

  dataCleanedJoinedWide_2_mixed <-
    full_join(
      dataCleanedJoinedWide_2_author,
      dataCleanedJoinedWide_2_date
    ) %>%
    full_join(dataCleanedJoinedWide_2_journal) %>%
    rowwise() %>%
    mutate(
      referenceTranslated_scoreComplement_total =
        sum(
          referenceTranslated_scoreComplement_date,
          referenceTranslated_scoreComplement_author,
          referenceTranslated_scoreComplement_journal,
          na.rm = TRUE
        )
    )

  dataTranslated_2_tmp <- dataTranslated_2 %>%
    select(
      organismType,
      organismValue,
      organismDetected,
      referenceType,
      referenceValue,
      referenceTranslated_title,
      referenceTranslated_scoreCrossref,
      level
    ) %>%
    left_join(dataTranslated_title_distinct)

  dataTranslated_2_tmp_organism <- dataTranslated_2_tmp %>%
    filter(!is.na(referenceTranslated_scoreTitleOrganism)) %>%
    left_join(dataCleanedJoinedWide_2_mixed) %>%
    group_by(
      organismType,
      organismValue,
      organismDetected,
      referenceValue
    ) %>%
    arrange(desc(as.numeric(referenceTranslated_scoreCrossref))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreComplement_total
    ))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreTitleOrganism
    ))) %>%
    ungroup() %>%
    distinct(organismType,
      organismValue,
      organismDetected,
      referenceValue,
      .keep_all = TRUE
    ) %>%
    select(-level)

  dataTranslated_2_tmp_no_organism <- dataTranslated_2_tmp %>%
    filter(is.na(referenceTranslated_scoreTitleOrganism)) %>%
    left_join(dataCleanedJoinedWide_2_mixed) %>%
    arrange(desc(as.numeric(referenceTranslated_scoreCrossref))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreComplement_total
    ))) %>%
    ungroup() %>%
    distinct(organismType,
      organismValue,
      organismDetected,
      referenceValue,
      .keep_all = TRUE
    ) %>%
    select(-level)

  dataTranslated_full <- bind_rows(
    dataCleanedJoinedWide_1 %>%
      left_join(dataTranslated_1),
    dataTranslated_2_tmp_no_organism %>%
      left_join(dataTranslated_2),
    dataTranslated_2_tmp_organism %>%
      left_join(dataTranslated_2),
    dataTranslated_3
  )

  dataCleanedJoinedWideScore <- dataTranslated %>%
    left_join(dataTranslated_full) %>%
    group_by(
      organismType,
      organismValue,
      organismDetected,
      referenceType,
      referenceValue
    ) %>%
    arrange(desc(as.numeric(referenceTranslated_scoreCrossref))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreComplement_total
    ))) %>%
    arrange(desc(as.numeric(
      referenceTranslated_scoreTitleOrganism
    ))) %>%
    ungroup() %>%
    distinct(organismType,
      organismValue,
      organismDetected,
      referenceValue,
      .keep_all = TRUE
    ) %>%
    mutate(referenceTranslated_doi = toupper(referenceTranslated_doi))

  subDataClean_doi <- dataCleanedJoinedWideScore %>%
    filter(!is.na(referenceTranslated_doi)) %>%
    distinct(referenceTranslated_doi) %>%
    mutate_all(as.character)

  subDataClean_pmid <- dataCleanedJoinedWideScore %>%
    filter(referenceType == "pubmed") %>%
    distinct(referenceValue) %>%
    mutate_all(as.character)

  # log_debug("adding PMID and PMCID")
  df_doi <- left_join(subDataClean_doi,
    PMC_ids,
    by = c("referenceTranslated_doi" = "DOI")
  ) %>%
    filter(!is.na(PMID) | !is.na(PMCID)) %>%
    select(
      referenceTranslated_doi,
      referenceTranslated_pmid = PMID,
      referenceTranslated_pmcid = PMCID
    )

  df_pubmed <- left_join(subDataClean_pmid,
    PMC_ids,
    by = c("referenceValue" = "PMID")
  ) %>%
    filter(!is.na(referenceValue) | !is.na(PMCID)) %>%
    select(referenceValue,
      referenceTranslated_pmcid = PMCID
    ) %>%
    mutate(referenceTranslated_pmid = referenceValue)

  tableJoined <- left_join(dataCleanedJoinedWideScore, df_doi)

  referenceTable <-
    left_join(tableJoined,
      df_pubmed,
      by = c("referenceValue" = "referenceValue")
    ) %>%
    mutate(
      referenceTranslated_pmid = ifelse(
        test = !is.na(referenceTranslated_pmid.x),
        yes = referenceTranslated_pmid.x,
        no = referenceTranslated_pmid.y
      ),
      referenceTranslated_pmcid = ifelse(
        test = !is.na(referenceTranslated_pmcid.x),
        yes = referenceTranslated_pmcid.x,
        no = referenceTranslated_pmcid.y
      )
    ) %>%
    select(
      organismType,
      organismValue,
      organismDetected,
      referenceType,
      referenceValue,
      referenceCleanedDoi = referenceTranslated_doi,
      referenceCleanedPmcid = referenceTranslated_pmcid,
      referenceCleanedPmid = referenceTranslated_pmid,
      referenceCleanedTitle = referenceTranslated_title,
      referenceCleaned_journal = referenceTranslated_journal,
      referenceCleaned_date = referenceTranslated_date,
      referenceCleaned_author = referenceTranslated_author,
      referenceCleaned_score_crossref = referenceTranslated_scoreCrossref,
      referenceCleaned_score_distance = referenceTranslated_scoreDistance,
      referenceCleaned_score_titleOrganism = referenceTranslated_scoreTitleOrganism,
      referenceCleaned_score_complementDate = referenceTranslated_scoreComplement_date,
      referenceCleaned_score_complementAuthor = referenceTranslated_scoreComplement_author,
      referenceCleaned_score_complementJournal = referenceTranslated_scoreComplement_journal,
      referenceCleaned_score_complementTotal = referenceTranslated_scoreComplement_total
    ) %>%
    distinct() %>%
    mutate(across(everything(), ~ y_as_na(.x, "NULL")))

  # log_debug("exporting ...")
  # log_debug(pathDataInterimTablesProcessedReferenceFile)
  write_delim(
    x = referenceTable,
    file = outpath,
    delim = "\t",
    na = ""
  )
}

log_debug("joining results together")
data3 <- do.call(
  "rbind",
  lapply(
    list.files(
      path = file.path(pathDataInterimTablesTranslatedReference),
      pattern = "*_processed.tsv.gz",
      full.names = FALSE
    ),
    function(x) {
      read_delim(
        file = gzfile(file.path(
          pathDataInterimTablesTranslatedReference, x
        )),
        delim = "\t",
        escape_double = TRUE,
        trim_ws = TRUE
      ) %>%
        mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedReferenceFile)
write_delim(
  x = data3,
  delim = "\t",
  file = pathDataInterimTablesProcessedReferenceFile,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
