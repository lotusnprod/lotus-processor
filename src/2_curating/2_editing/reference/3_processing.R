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
  readr::read_delim(
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
    ) |>
    dplyr::mutate_all(as.character) |>
    dplyr::tibble()
} else if (mode == "custom") {
  PMC_ids <-
    data.frame(
      DOI = NA,
      PMCID = NA,
      PMID = NA
    ) |>
    dplyr::mutate_all(as.character) |>
    dplyr::tibble()
} else {
  log_debug("loading pmcid file, this may take a while")
  PMC_ids <- readr::read_delim(
    file = pathDataExternalTranslationSourcePubmedFile,
    delim = ",",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
    dplyr::filter(!is.na(DOI) | !is.na(PMID)) |>
    dplyr::select(
      DOI,
      PMCID,
      PMID
    ) |>
    dplyr::mutate(DOI = toupper(DOI)) |>
    dplyr::mutate_all(as.character) |>
    dplyr::tibble()
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
create_dir(export = pathDataInterimTablesProcessedReference)

for (i in num) {
  inpath <-
    paste0(
      pathDataInterimTablesTranslatedReference,
      "/",
      stringr::str_pad(
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
      stringr::str_pad(
        string = i,
        width = 6,
        pad = "0"
      ),
      "_processed.tsv.gz"
    )

  log_debug(paste("step", i / cut, "of", length))

  # log_debug("loading translated file")
  dataTranslated <-
    readr::read_delim(
      file = inpath,
      delim = "\t",
      col_types = cols(.default = "c"),
      locale = locales
    )

  dataCrossref_step <- dataCrossref |>
    dplyr::semi_join(
      dataTranslated,
      by = c(
        "referenceOriginal" = "referenceValue",
        "origin" = "referenceType"
      )
    ) |>
    dplyr::distinct(
      referenceOriginal,
      referenceTranslatedType,
      origin,
      level,
      .keep_all = TRUE
    ) |>
    dplyr::semi_join(
      dataTranslated,
      by = c(
        "referenceOriginal" = "referenceValue",
        "origin" = "referenceType"
      )
    ) |>
    dplyr::distinct(
      referenceOriginal,
      referenceTranslatedType,
      origin,
      level,
      .keep_all = TRUE
    )
  if (nrow(dataCrossref_step) != 0) {
    # log_debug("selecting minimal crossref data for DOI, PUBMED, TITLE")
    dataCrossref_1 <- dataCrossref_step |>
      dplyr::filter(
        referenceTranslatedType == "doi" |
          referenceTranslatedType == "title" |
          referenceTranslatedType == "scoreCrossref" |
          referenceTranslatedType == "scoreDistance"
      ) |>
      tidyr::pivot_wider(
        names_from = referenceTranslatedType,
        names_prefix = "referenceTranslated_",
        values_from = referenceTranslatedValue
      )

    # log_debug("selecting minimal crossref data for ORIGINAL, PUBDETAILS, SPLIT")
    dataCrossref_2 <- dataCrossref_step |>
      dplyr::filter(
        referenceTranslatedType == "doi" |
          referenceTranslatedType == "title" |
          referenceTranslatedType == "scoreCrossref" |
          referenceTranslatedType == "journal" |
          referenceTranslatedType == "date" |
          referenceTranslatedType == "author"
      ) |>
      tidyr::pivot_wider(
        names_from = referenceTranslatedType,
        names_prefix = "referenceTranslated_",
        values_from = referenceTranslatedValue
      )

    # log_debug("selecting minimal data for organism in title recognition")
    dataTranslated_title <- dataTranslated |>
      dplyr::filter(!is.na(organismDetected)) |>
      dplyr::inner_join(
        dataCrossref_step |>
          dplyr::filter(referenceTranslatedType == "title"),
        by = c(
          "referenceValue" = "referenceOriginal",
          "referenceType" = "origin"
        )
      )

    dataTranslated_title_distinct <- dataTranslated_title |>
      dplyr::distinct(organismDetected, referenceTranslatedValue)

    # log_debug("checking for organism in title, may take a while if running full mode")
    dataTranslated_title_distinct <-
      dataTranslated_title_distinct |>
      dplyr::mutate(
        referenceTranslated_scoreTitleOrganism = ifelse(
          test = stringr::str_detect(
            string = fixed(tolower(referenceTranslatedValue)),
            pattern = fixed(tolower(organismDetected))
          ),
          yes = 1,
          no = 0
        )
      ) |>
      dplyr::filter(referenceTranslated_scoreTitleOrganism == 1) |>
      dplyr::distinct(
        organismDetected,
        referenceTranslated_title = referenceTranslatedValue,
        referenceTranslated_scoreTitleOrganism
      )

    # log_debug("splitting cases")
    dataTranslated_1 <- dataTranslated |>
      dplyr::filter(
        referenceType == "doi" |
          referenceType == "pubmed" |
          referenceType == "title"
      ) |>
      dplyr::left_join(
        dataCrossref_1,
        by = c(
          "referenceValue" = "referenceOriginal",
          "referenceType" = "origin"
        )
      )

    dataTranslated_2 <- dataTranslated |>
      dplyr::filter(
        referenceType == "original" |
          referenceType == "publishingDetails" |
          referenceType == "split"
      ) |>
      dplyr::left_join(
        dataCrossref_2,
        by = c(
          "referenceValue" = "referenceOriginal",
          "referenceType" = "origin"
        )
      )

    dataTranslated_3 <- dataTranslated |>
      dplyr::filter(
        referenceType != "doi" &
          referenceType != "pubmed" &
          referenceType != "title" &
          referenceType != "original" &
          referenceType != "publishingDetails" &
          referenceType != "split"
      )

    dataCleanedJoinedWide_1 <- dataTranslated_1 |>
      dplyr::distinct(
        organismType,
        organismValue,
        organismDetected,
        referenceType,
        referenceValue,
        referenceTranslated_title,
        referenceTranslated_scoreCrossref,
        referenceTranslated_scoreDistance,
        level
      ) |>
      dplyr::left_join(dataTranslated_title_distinct) |>
      dplyr::select(-referenceTranslated_title) |>
      dplyr::distinct() |>
      dplyr::group_by(
        organismType,
        organismValue,
        organismDetected,
        referenceType,
        referenceValue
      ) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreCrossref
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreTitleOrganism
      ))) |>
      dplyr::arrange(as.numeric(referenceTranslated_scoreDistance)) |>
      dplyr::ungroup() |>
      dplyr::distinct(
        organismType,
        organismValue,
        organismDetected,
        referenceValue,
        .keep_all = TRUE
      ) |>
      dplyr::select(-level)

    dataCleanedJoinedWide_2_author <- dataTranslated_2 |>
      dplyr::distinct(
        referenceType,
        referenceValue,
        referenceTranslated_doi,
        referenceTranslated_title,
        referenceTranslated_author
      ) |>
      dplyr::mutate(
        referenceTranslated_scoreComplement_author = ifelse(
          stringr::str_detect(
            string = tolower(referenceValue),
            pattern = fixed(tolower(word(
              referenceTranslated_author,
              1
            )))
          ),
          yes = 1,
          no = 0
        )
      ) |>
      dplyr::filter(referenceTranslated_scoreComplement_author == 1)

    dataCleanedJoinedWide_2_date <- dataTranslated_2 |>
      dplyr::distinct(
        referenceType,
        referenceValue,
        referenceTranslated_doi,
        referenceTranslated_title,
        referenceTranslated_date
      ) |>
      dplyr::mutate(
        referenceTranslated_year = substr(referenceTranslated_date, 1, 4),
        referenceTranslated_scoreComplement_date = dplyr::if_else(
          !is.na(referenceTranslated_year) &
            stringr::str_detect(referenceValue, referenceTranslated_year),
          1L,
          0L
        )
      ) |>
      dplyr::filter(referenceTranslated_scoreComplement_date == 1)

    dataCleanedJoinedWide_2_journal <- dataTranslated_2 |>
      dplyr::distinct(
        referenceType,
        referenceValue,
        referenceTranslated_doi,
        referenceTranslated_title,
        referenceTranslated_journal
      ) |>
      dplyr::mutate(
        referenceTranslated_scoreComplement_journal = ifelse(
          stringr::str_detect(
            string = tolower(referenceValue),
            pattern = stringr::fixed(tolower(referenceTranslated_journal))
          ),
          yes = 1,
          no = 0
        )
      ) |>
      dplyr::filter(referenceTranslated_scoreComplement_journal == 1)

    dataCleanedJoinedWide_2_mixed <-
      dplyr::full_join(
        dataCleanedJoinedWide_2_author,
        dataCleanedJoinedWide_2_date
      ) |>
      dplyr::full_join(dataCleanedJoinedWide_2_journal) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        referenceTranslated_scoreComplement_total = sum(
          referenceTranslated_scoreComplement_date,
          referenceTranslated_scoreComplement_author,
          referenceTranslated_scoreComplement_journal,
          na.rm = TRUE
        )
      )

    dataTranslated_2_tmp <- dataTranslated_2 |>
      dplyr::select(
        organismType,
        organismValue,
        organismDetected,
        referenceType,
        referenceValue,
        referenceTranslated_title,
        referenceTranslated_scoreCrossref,
        level
      ) |>
      dplyr::left_join(dataTranslated_title_distinct)

    dataTranslated_2_tmp_organism <- dataTranslated_2_tmp |>
      dplyr::filter(!is.na(referenceTranslated_scoreTitleOrganism)) |>
      dplyr::left_join(dataCleanedJoinedWide_2_mixed) |>
      dplyr::group_by(
        organismType,
        organismValue,
        organismDetected,
        referenceValue
      ) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreCrossref
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreComplement_total
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreTitleOrganism
      ))) |>
      dplyr::ungroup() |>
      dplyr::distinct(
        organismType,
        organismValue,
        organismDetected,
        referenceValue,
        .keep_all = TRUE
      ) |>
      dplyr::select(-level)

    dataTranslated_2_tmp_no_organism <- dataTranslated_2_tmp |>
      dplyr::filter(is.na(referenceTranslated_scoreTitleOrganism)) |>
      dplyr::left_join(dataCleanedJoinedWide_2_mixed) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreCrossref
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreComplement_total
      ))) |>
      dplyr::ungroup() |>
      dplyr::distinct(
        organismType,
        organismValue,
        organismDetected,
        referenceValue,
        .keep_all = TRUE
      ) |>
      dplyr::select(-level)

    dataTranslated_full <- dplyr::bind_rows(
      dataCleanedJoinedWide_1 |>
        dplyr::left_join(dataTranslated_1),
      dataTranslated_2_tmp_no_organism |>
        dplyr::left_join(dataTranslated_2),
      dataTranslated_2_tmp_organism |>
        dplyr::left_join(dataTranslated_2),
      dataTranslated_3
    )

    dataCleanedJoinedWideScore <- dataTranslated |>
      dplyr::left_join(dataTranslated_full) |>
      dplyr::group_by(
        organismType,
        organismValue,
        organismDetected,
        referenceType,
        referenceValue
      ) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreCrossref
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreComplement_total
      ))) |>
      dplyr::arrange(dplyr::desc(as.numeric(
        referenceTranslated_scoreTitleOrganism
      ))) |>
      dplyr::ungroup() |>
      dplyr::distinct(
        organismType,
        organismValue,
        organismDetected,
        referenceValue,
        .keep_all = TRUE
      ) |>
      dplyr::mutate(referenceTranslated_doi = toupper(referenceTranslated_doi))

    subDataClean_doi <- dataCleanedJoinedWideScore |>
      dplyr::filter(!is.na(referenceTranslated_doi)) |>
      dplyr::distinct(referenceTranslated_doi) |>
      dplyr::mutate_all(as.character)

    subDataClean_pmid <- dataCleanedJoinedWideScore |>
      dplyr::filter(referenceType == "pubmed") |>
      dplyr::distinct(referenceValue) |>
      dplyr::mutate_all(as.character)

    # log_debug("adding PMID and PMCID")
    df_doi <- dplyr::left_join(
      subDataClean_doi,
      PMC_ids,
      by = c("referenceTranslated_doi" = "DOI")
    ) |>
      dplyr::filter(!is.na(PMID) | !is.na(PMCID)) |>
      dplyr::select(
        referenceTranslated_doi,
        referenceTranslated_pmid = PMID,
        referenceTranslated_pmcid = PMCID
      )

    df_pubmed <- dplyr::left_join(
      subDataClean_pmid,
      PMC_ids,
      by = c("referenceValue" = "PMID")
    ) |>
      dplyr::filter(!is.na(referenceValue) | !is.na(PMCID)) |>
      dplyr::select(referenceValue, referenceTranslated_pmcid = PMCID) |>
      dplyr::mutate(referenceTranslated_pmid = referenceValue)

    tableJoined <-
      dplyr::left_join(dataCleanedJoinedWideScore, df_doi)

    referenceTable <-
      dplyr::left_join(
        tableJoined,
        df_pubmed,
        by = c("referenceValue" = "referenceValue")
      ) |>
      dplyr::mutate(
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
      ) |>
      dplyr::select(
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
      ) |>
      dplyr::distinct() |>
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ y_as_na(.x, "NULL")))

    # log_debug("exporting ...")
    # log_debug(pathDataInterimTablesProcessedReferenceFile)
    readr::write_delim(
      x = referenceTable,
      file = outpath,
      delim = "\t",
      na = ""
    )
  }
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
      readr::read_delim(
        file = gzfile(
          description = file.path(
            pathDataInterimTablesTranslatedReference,
            x
          )
        ),
        delim = "\t",
        escape_double = TRUE,
        trim_ws = TRUE
      ) |>
        dplyr::mutate_all(as.character)
    }
  )
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesProcessedReferenceFile)
readr::write_delim(
  x = data3,
  delim = "\t",
  file = pathDataInterimTablesProcessedReferenceFile,
  na = ""
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
