## this is very dirty for now, I'll be cleaning it later on
source("r/log_debug.R")
log_debug(
  "This script aims to establish filtering criteria to validate",
  "documented pairs according to manually analyzed ones"
)

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(plotly)
library(readr)
library(stringi)
library(tidyr)

log_debug("... functions")
source("r/filter_dirty.R")
source("r/myDirtyValidationFig.R")

log_debug("loading files ...")
oldDbNames <-
  readr::read_delim(
    file = "../data/external/dictionarySource/dbNames.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_AR_old <-
  readr::read_delim(
    file = "../data/validation/old/AR.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "AR")

sampleAllONPDB_PMA_old <-
  readr::read_delim(
    file = "../data/validation/old/PMA.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "PMA")

sampleAllONPDB_AR <-
  readr::read_delim(
    file = "../data/validation/new/AR.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "AR")

sampleAllONPDB_JB <-
  readr::read_delim(
    file = "../data/validation/new/JB.csv",
    delim = ",",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "JB") |>
  dplyr::select(1:20)

sampleAllONPDB_PMA <-
  readr::read_delim(
    file = "../data/validation/new/PM.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "PMA")

sampleAllONPDB_publishingDetails <-
  readr::read_delim(
    file = "../data/validation/new/publishingDetails.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_additionalSet <-
  readr::read_delim(
    file = "../data/validation/new/additionalSet.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleAllONPDB_additionalSetBis <-
  readr::read_delim(
    file = "../data/validation/new/additionalSetBis.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  )

sampleCondifent_PMA <-
  readr::read_delim(
    file = "../data/validation/confident/100confidentPMAChecked.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(curator == "PMA") |>
  dplyr::mutate(curator = "PMA2")

log_debug("... documented pairs")
inhouseDbMinimal <-
  readr::read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("... structure metadata")
structureMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesStructureMetadata,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct(
    structureCleanedSmiles,
    structureCleaned_smiles2D,
    structureCleanedInchi,
    structureCleaned_inchi2D,
    structureCleanedInchikey,
    structureCleaned_inchikey2D,
    structureCleaned_molecularFormula,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    structureCleaned_cid,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    structureCleaned_nameIupac = dplyr::if_else(
      condition = is.na(structureCleaned_nameIupac),
      true = structureCleanedInchikey,
      false = structureCleaned_nameIupac
    ),
    structureCleaned_nameTraditional = dplyr::if_else(
      condition = is.na(structureCleaned_nameTraditional),
      true = structureCleanedInchikey,
      false = structureCleaned_nameTraditional
    )
  )

log_debug("... organism metadata")
organismMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesOrganismMetadata,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct(
    organismCleaned,
    organismCleaned_id,
    organismCleaned_dbTaxo,
    # organismCleaned_dbTaxoTaxonIds,
    organismCleaned_dbTaxoTaxonRanks,
    organismCleaned_dbTaxoTaxonomy
  )

log_debug("... reference metadata")
referenceMetadata <-
  readr::read_delim(
    file = pathDataInterimDictionariesReferenceMetadata,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  ) |>
  dplyr::distinct(
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

table_old <- sampleAllONPDB_old |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate(
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
  ) |>
  tidyr::pivot_longer(
    6:10,
    names_to = c("origin1", "referenceType"),
    names_sep = "_",
    values_to = "referenceValue",
    values_drop_na = TRUE
  ) |>
  tidyr::pivot_longer(
    3:5,
    names_to = c("origin2", "structureType"),
    names_sep = "_",
    values_to = "structureValue",
    values_drop_na = TRUE
  ) |>
  dplyr::filter(!is.na(referenceValue)) |>
  dplyr::group_by(inchikeySanitized) |>
  dplyr::distinct(referenceCleanedDoi, .keep_all = TRUE) |>
  dplyr::ungroup() |>
  dplyr::select(-origin1, -origin2) |>
  dplyr::select(
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

sampleAllONPDB <- dplyr::bind_rows(
  sampleAllONPDB_AR,
  sampleAllONPDB_JB,
  sampleAllONPDB_PMA,
  sampleAllONPDB_publishingDetails,
  sampleAllONPDB_additionalSet,
  sampleAllONPDB_additionalSetBis,
  sampleCondifent_PMA
)

table <- sampleAllONPDB |>
  dplyr::mutate_all(as.character) |>
  dplyr::mutate(
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
  ) |>
  dplyr::select(
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

globalSample <- dplyr::bind_rows(table_old, table) |>
  dplyr::mutate(organismValue = organismOriginal) |>
  dplyr::distinct(
    database,
    # organismValue,
    # structureType,
    # structureValue,
    # referenceType,
    # referenceValue,
    organismCleaned,
    structureCleanedInchi,
    structureCleanedSmiles,
    structureCleanedInchikey,
    referenceCleanedDoi,
    curator,
    validated,
    comments
  ) |>
  dplyr::filter(!is.na(validated)) |>
  dplyr::mutate(referenceCleanedDoi = toupper(referenceCleanedDoi)) |>
  dplyr::left_join(inhouseDbMinimal) |>
  dplyr::select(
    -database,
    -referenceCleanedTitle,
    -organismType,
    -organismValue,
    -structureType,
    -structureValue,
    -referenceType,
    -referenceValue,
  ) |>
  dplyr::distinct()

a <- paste0("\\b", oldDbNames$oldDbName, "\\b")
b <- oldDbNames$newDbName

log_debug("adding metadata")
inhouseDbFull <- inhouseDbMinimal |>
  dplyr::left_join(referenceMetadata)

log_debug("joining manual validation results with documented pairs")
realMetaSample <- dplyr::inner_join(globalSample, inhouseDbFull) |>
  dplyr::distinct(
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
  ) |>
  dplyr::filter(
    referenceType %in%
      c(
        "doi",
        "pubmed",
        "title",
        "publishingDetails",
        "split",
        "original"
      )
  )

if (mode == "full") {
  log_debug("filtering results ...")
  log_debug("... validated set")
  realSampleFiltered <-
    filter_dirty(dataframe = realMetaSample)

  log_debug("... rejected set")
  antiFilter <- dplyr::anti_join(realMetaSample, realSampleFiltered)

  log_debug("counting results ...")
  log_debug("... per category on validation set")
  table_count <- myDirtyF(table = realMetaSample)

  log_debug("... global on validation set")
  table_count_global <- myDirtyC(table = realMetaSample)

  log_debug("... per category on validated set")
  tableFiltered_count <-
    myDirtyF(table = realSampleFiltered)

  log_debug("... global on validated set")
  tableFiltered_count_global <-
    myDirtyC(table = realSampleFiltered)

  log_debug("... per category on rejected set")
  tableAntiFiltered_count <-
    myDirtyF(table = antiFilter)

  log_debug("... global on rejected set")
  tableAntiFiltered_count_global <-
    myDirtyC(table = antiFilter)

  log_debug("visualizing ...")
  log_debug("... validation set per category")
  fig_full <-
    myDirtyP(
      table = table_count,
      yaxismax = 300,
      title = "full version"
    )
  fig_full

  log_debug("... validated set per category")
  fig_filtered <-
    myDirtyP(
      table = tableFiltered_count,
      yaxismax = 300,
      title = "filtered version"
    )
  fig_filtered

  log_debug("... rejected set per category")
  fig_anti <-
    myDirtyP(
      table = tableAntiFiltered_count,
      yaxismax = 140,
      title = "anti version"
    )
  fig_anti

  log_debug("... validation set global")
  newfull <- myDirtyQ(
    table = table_count_global,
    yaxismax = 550,
    title = "new full version"
  )
  newfull

  log_debug("... validated set global")
  newfiltered <- myDirtyQ(
    table = tableFiltered_count_global,
    yaxismax = 550,
    title = "new filtered version"
  )
  newfiltered

  log_debug("... rejected set global")
  antifull <- myDirtyQ(
    table = tableAntiFiltered_count_global,
    yaxismax = 100,
    title = "anti full version"
  )
  antifull

  log_debug("calculating statistics ...")
  old <- table_count |>
    dplyr::select(referenceType, tot, y, n, ratio1 = ratio)

  log_debug("... true positives and false positives")
  new <- tableFiltered_count |>
    dplyr::select(
      referenceType,
      fil = tot,
      tp = y,
      fp = n,
      ratio2 = ratio
    )

  log_debug("... true negatives and false negatives")
  anti <- tableAntiFiltered_count |>
    dplyr::select(
      referenceType,
      anti = tot,
      fn = y,
      tn = n,
      ratioAnti = ratio
    )

  f1Table <- dplyr::full_join(old, new) %>%
    replace(is.na(.), 0)

  beta <- 0.5

  log_debug("... precision, recall and Fbeta", beta, "score")
  f1Table <- dplyr::full_join(f1Table, anti) %>%
    dplyr::mutate(
      tpfn = tp + fn,
      tpfp = tp + fp
    ) %>%
    dplyr::mutate(
      recall = tp / tpfn,
      precision = tp / tpfp
    ) %>%
    dplyr::mutate(
      rxp = recall * precision,
      rpp = recall + precision
    ) %>%
    dplyr::mutate(f2 = 2 * rxp / rpp) %>%
    dplyr::select(-tpfn, -tpfp, -rxp, -rpp) %>%
    dplyr::mutate(
      fbeta = (1 + beta^2) *
        (precision * recall) /
        ((precision * beta^2) + recall)
    ) %>%
    replace(. == "NaN", 0)

  "%ni%" <- Negate("%in%")

  log_debug("... correcting with references ratios")
  refRatios <- inhouseDbFull %>%
    dplyr::filter(database %ni% forbidden_export) %>%
    dplyr::group_by(referenceType) %>%
    dplyr::count() %>%
    dplyr::mutate(prop = n / sum(.$n, na.rm = TRUE)) %>%
    dplyr::select(-n)

  f2Table <- dplyr::left_join(f1Table, refRatios) %>%
    dplyr::mutate(
      f2corr_temp = f2 * prop,
      recall_temp = recall * prop,
      precision_temp = precision * prop
    ) %>%
    dplyr::mutate(
      f2corr = sum(.$f2corr_temp, na.rm = TRUE),
      recallcorr = sum(.$recall_temp, na.rm = TRUE),
      precisioncorr = sum(.$precision_temp, na.rm = TRUE)
    ) %>%
    dplyr::mutate()

  f2Table_reworked <- f2Table %>%
    dplyr::arrange(dplyr::desc(prop)) %>%
    dplyr::select(
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
      dplyr::tibble(
        `reference type` = "Total",
        `true positives` = sum(.$`true positives`, na.rm = TRUE),
        `false positives` = sum(.$`false positives`, na.rm = TRUE),
        `false negatives` = sum(.$`false negatives`, na.rm = TRUE),
        `true negatives` = sum(.$`true negatives`, na.rm = TRUE),
        `relative abundance` = sum(.$`relative abundance`, na.rm = TRUE),
      ),
      dplyr::tibble(
        `reference type` = "Corrected total",
        `precision` = unique(.$precisioncorr),
        `recall` = unique(.$recallcorr),
        `F0.5 score` = unique(.$f2corr),
      )
    ) %>%
    dplyr::select(
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
}

log_debug(
  "applying the filtering criteria to the whole DB, this may take a while"
)
openDb <- inhouseDbFull |>
  filter_dirty() |>
  dplyr::left_join(structureMetadata) |>
  dplyr::left_join(organismMetadata) |>
  dplyr::distinct(
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
  ) |>
  dplyr::mutate(
    structureCleaned_inchikey2D = substring(
      text = structureCleanedInchikey,
      first = 1,
      last = 14
    )
  ) |>
  dplyr::select(
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
    structureCleaned_molecularFormula,
    structureCleaned_cid,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

log_debug("outputting closed pairs")
closedDb <- inhouseDbFull |>
  dplyr::filter(database %in% forbidden_export) |>
  dplyr::left_join(structureMetadata) |>
  dplyr::left_join(organismMetadata) |>
  dplyr::distinct(
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
    structureCleaned_cid,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle,
    .keep_all = TRUE
  ) |>
  mutate(
    structureCleaned_inchikey2D = substring(
      text = structureCleanedInchikey,
      first = 1,
      last = 14
    )
  ) |>
  dplyr::select(
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
    structureCleaned_molecularFormula,
    structureCleaned_cid,
    structureCleaned_nameIupac,
    structureCleaned_nameTraditional,
    structureCleaned_stereocenters_total,
    structureCleaned_stereocenters_unspecified,
    referenceCleanedDoi,
    referenceCleanedPmcid,
    referenceCleanedPmid,
    referenceCleanedTitle
  )

log_debug("outputing correct entries from manually validated set")
manuallyValidatedSet <- realMetaSample |>
  dplyr::filter(validated == "Y") |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedTitle
  )

log_debug("outputing incorrect entries from validated set")
manuallyRemovedEntries <- realMetaSample |>
  dplyr::filter(validated != "Y") |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedTitle
  )

openDbClean <- dplyr::anti_join(openDb, manuallyRemovedEntries)

if (mode == "full") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  validationSet <- dplyr::anti_join(openDbClean, realMetaSample) |>
    dplyr::sample_n(100)
}

log_debug("loading validation set")
validationSetFilled_1 <-
  readr::read_delim(
    file = "../data/validation/validationSet.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(!is.na(validated)) |>
  dplyr::mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

if (mode == "full") {
  set.seed(
    seed = 42,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  validationSet2 <-
    dplyr::anti_join(openDbClean, validationSetFilled_1) |>
    dplyr::sample_n(13)
}

log_debug("loading validation set bis")
validationSetFilled_2 <-
  readr::read_delim(
    file = "../data/validation/validationSetBis.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(!is.na(validated)) |>
  dplyr::mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

if (mode == "test") {
  validationSet3 <-
    dplyr::anti_join(openDbClean, validationSetFilled_1) |>
    dplyr::anti_join(validationSetFilled_2)
}

log_debug("loading validation set ter")
validationSetFilled_3 <-
  readr::read_delim(
    file = "../data/validation/validationSetTer.tsv",
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales,
    escape_double = FALSE,
    trim_ws = TRUE
  ) |>
  dplyr::filter(!is.na(validated)) |>
  dplyr::mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

# log_debug("loading validation set tetra")
# validationSetFilled_4 <-
#   readr::read_delim(
#     file = "../data/validation/validationSetTetra.tsv",
#     delim = "\t",
#     col_types = cols(.default = "c"),
#     locale = locales,
#     escape_double = FALSE,
#     trim_ws = TRUE
#   ) |>
#   dplyr::filter(!is.na(validated)) |>
#   dplyr::mutate(referenceCleanedDoi = toupper(referenceCleanedDoi))

## add validationSetFilled_4 below in case

validationSetFilled <-
  dplyr::bind_rows(
    validationSetFilled_1,
    validationSetFilled_2,
    validationSetFilled_3
  ) |>
  dplyr::mutate(organismValue = organismOriginal) |>
  dplyr::select(
    structureCleanedInchikey = structureCleanedInchikey3D,
    dplyr::everything()
  ) |>
  dplyr::distinct(
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
  ) |>
  dplyr::left_join(inhouseDbMinimal) |>
  dplyr::select(-referenceCleanedTitle, -organismType) |>
  dplyr::distinct()

validationSetFilled$database <- stringi::stri_replace_all_regex(
  str = validationSetFilled$database,
  pattern = a,
  replacement = b,
  case_insensitive = FALSE,
  vectorize_all = FALSE
)

realValidationSetFilled <-
  dplyr::inner_join(validationSetFilled, openDbClean) |>
  dplyr::distinct(
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

finalStats <- realValidationSetFilled |>
  dplyr::group_by(referenceType) |>
  dplyr::count(validated == "Y") |>
  dplyr::ungroup()

finalStats_reworked <- finalStats %>%
  dplyr::mutate(
    `true positives` = dplyr::if_else(
      condition = `validated == "Y"` == TRUE,
      true = n,
      false = 0
    ),
    `false positives` = dplyr::if_else(
      condition = `validated == "Y"` == FALSE,
      true = n,
      false = 0
    )
  ) %>%
  dplyr::select(
    `reference type` = referenceType,
    dplyr::everything(),
    -n,
    -`validated == "Y"`
  ) %>%
  dplyr::group_by(`reference type`) %>%
  dplyr::summarize(
    `true positives Validation` = max(`true positives`),
    `false positives Validation` = max(`false positives`)
  ) %>%
  rbind(
    .,
    tibble(
      `reference type` = "Total",
      `true positives Validation` = sum(
        .$`true positives Validation`,
        na.rm = TRUE
      ),
      `false positives Validation` = sum(
        .$`false positives Validation`,
        na.rm = TRUE
      )
    )
  )

if (mode == "full") {
  finalTable <-
    dplyr::left_join(f2Table_reworked, finalStats_reworked) %>%
    dplyr::mutate_if(.predicate = is.numeric, ~ round(., digits = 2))
}

log_debug("outputing correct entries from manually validated set")
manuallyValidatedSet2 <- realValidationSetFilled |>
  dplyr::filter(validated == "Y") |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedTitle
  )

manuallyValidatedSet3 <-
  dplyr::bind_rows(manuallyValidatedSet, manuallyValidatedSet2)

log_debug("outputing incorrect entries from validated set")
manuallyRemovedEntries2 <- realValidationSetFilled |>
  dplyr::filter(validated != "Y")

manuallyRemovedEntries3 <-
  dplyr::bind_rows(manuallyRemovedEntries, manuallyRemovedEntries2) |>
  dplyr::distinct(
    organismCleaned,
    structureCleanedInchi,
    structureCleanedInchikey,
    structureCleanedSmiles,
    referenceCleanedDoi,
    referenceCleanedTitle
  )

openDbClean2 <-
  dplyr::anti_join(openDbClean, manuallyRemovedEntries3) |>
  dplyr::filter(!database %in% forbidden_export)

log_debug("removing dimers")
## at the moment no solution for it so discarding for safety
openDbClean3 <- openDbClean2 |>
  dplyr::filter(
    !grepl(
      pattern = "\\.",
      x = structureCleanedSmiles
    )
  )

closedDb3 <- closedDb |>
  dplyr::filter(
    !grepl(
      pattern = "\\.",
      x = structureCleanedSmiles
    )
  )

log_debug("non-validated entries")
non_validated <- inhouseDbMinimal |>
  dplyr::anti_join(openDbClean3) |>
  dplyr::filter(!database %in% forbidden_export) |>
  dplyr::distinct(
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
    referenceCleanedTitle,
  )

log_debug("exporting")
create_dir(export = pathDataInterimTablesAnalyzed)

if (mode == "full") {
  log_debug("../data/validation/manuallyValidated.tsv.gz")

  write.table(
    x = manuallyValidatedSet3,
    file = gzfile(
      description = "../data/validation/manuallyValidated.tsv.gz",
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )

  log_debug("../data/validation/manuallyRemoved.tsv.gz")

  write.table(
    x = manuallyRemovedEntries3,
    file = gzfile(
      description = "../data/validation/manuallyRemoved.tsv.gz",
      compression = 9,
      encoding = "UTF-8"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )

  write.table(
    x = finalTable,
    file = "../data/validation/tableStats.csv",
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = ",",
    fileEncoding = "UTF-8"
  )
}

log_debug(
  nrow(
    openDbClean3 |>
      dplyr::distinct(
        structureCleanedInchikey,
        organismCleaned,
        referenceCleanedDoi
      )
  ),
  "referenced pairs are being exported to",
  pathDataInterimTablesAnalyzedPlatinum
)
readr::write_delim(
  x = openDbClean3,
  delim = "\t",
  file = pathDataInterimTablesAnalyzedPlatinum
)

log_debug(file.path(pathDataInterimTablesAnalyzed, "closed.tsv.gz"))
readr::write_delim(
  x = closedDb3,
  delim = "\t",
  file = file.path(pathDataInterimTablesAnalyzed, "closed.tsv.gz")
)

readr::write_delim(
  x = non_validated,
  delim = "\t",
  file = pathDataInterimTablesAnalyzedGarbage
)

if (exists("validationSet")) {
  log_debug(file.path(
    pathDataInterimTablesAnalyzed,
    "validationSet.tsv"
  ))
  write.table(
    x = validationSet,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "validationSet.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("validationSet2")) {
  log_debug(file.path(
    pathDataInterimTablesAnalyzed,
    "validationSetBis.tsv"
  ))

  write.table(
    x = validationSet2,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "validationSetBis.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("validationSet3")) {
  log_debug(file.path(
    pathDataInterimTablesAnalyzed,
    "validationSetTer.tsv"
  ))
  write.table(
    x = validationSet3,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "validationSetTer.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("validationSet4")) {
  log_debug(file.path(
    pathDataInterimTablesAnalyzed,
    "validationSetTetra.tsv"
  ))
  write.table(
    x = validationSet3,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "validationSetTetra.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
