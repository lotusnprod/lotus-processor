source("r/log_debug.R")
log_debug("This script samples some entries to then manually check their validity.")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

log_debug("... libraries")
library(dplyr)
library(readr)

log_debug("loading db, if running fullmode, this may take a while")
openDbMinimal <-
  readr::read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c"),
    locale = locales
  )

log_debug("sampling ...")
log_debug("... DOI")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "doi")) >= 30) {
  sampleONPDB_doi <- openDbMinimal |>
    dplyr::filter(referenceType == "doi") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = NA,
      validated = NA,
      comments = NA
    )
}

log_debug("... original")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "original")) >= 30) {
  sampleONPDB_original <- openDbMinimal |>
    dplyr::filter(referenceType == "original") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = NA,
      validated = NA,
      comments = NA
    )
}

log_debug("... PMID")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "pubmed")) >= 30) {
  sampleONPDB_pubmed <- openDbMinimal |>
    dplyr::filter(referenceType == "pubmed") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = NA,
      validated = NA,
      comments = NA
    )
}

log_debug("... split")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "split")) >= 30) {
  sampleONPDB_split <- openDbMinimal |>
    dplyr::filter(referenceType == "split") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = NA,
      validated = NA,
      comments = NA
    )
}

log_debug("... title")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "title")) >= 30) {
  sampleONPDB_title <- openDbMinimal |>
    dplyr::filter(referenceType == "title") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = NA,
      validated = NA,
      comments = NA
    )
}

log_debug("... publishing details")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "publishingDetails")) >= 30) {
  sampleONPDB_publishingDetails <- openDbMinimal |>
    dplyr::filter(referenceType == "publishingDetails") |>
    dplyr::sample_n(30) |>
    dplyr::mutate(
      curator = "AR",
      validated = NA,
      comments = NA
    )
}

# get0 if to avoid error in minimal mode if df not present
sampleONPDB <- dplyr::bind_rows(
  get0(x = "sampleONPDB_doi"),
  get0(x = "sampleONPDB_original"),
  get0(x = "sampleONPDB_publishingDetails"),
  get0(x = "sampleONPDB_pubmed"),
  get0(x = "sampleONPDB_split"),
  get0(x = "sampleONPDB_title")
) |>
  dplyr::mutate_all(as.character)

log_debug("... attributing curator")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
sampleONPDB <- sampleONPDB[sample(nrow(sampleONPDB)), ]

sampleONPDB[1:50, "curator"] <- "AR"

sampleONPDB[51:100, "curator"] <- "JB"

sampleONPDB[101:150, "curator"] <- "PMA"

log_debug("... knapsack entries")
# set.seed(seed = 42,
#          kind = "Mersenne-Twister",
#          normal.kind = "Inversion")
# sampleKnapsack <- openDbMinimal %>%
#   dplyr::filter(database == "kna_1") %>%
#   dplyr::sample_n(150) %>%
#   dplyr::mutate(
#     curator = sample(c("AR", "JB", "PMA"),
#                      size = nrow(.),
#                      replace = TRUE),
#     validated = NA,
#     comments = NA
#   )

log_debug("... additional entries")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal |>
  dplyr::filter(referenceType == "title")) >= 85 &
  nrow(openDbMinimal |>
    dplyr::filter(referenceType == "publishingDetails")) >= 25 &
  nrow(openDbMinimal |>
    dplyr::filter(referenceType == "split")) >= 41 &
  nrow(openDbMinimal |>
    dplyr::filter(referenceType == "publishingDetails")) >= 58) {
  additionalSet <-
    dplyr::bind_rows(
      A <- openDbMinimal |>
        dplyr::filter(referenceType == "title") |>
        dplyr::sample_n(85) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbMinimal |>
        dplyr::filter(referenceType == "publishingDetails") |>
        dplyr::sample_n(25) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbMinimal |>
        dplyr::filter(referenceType == "split") |>
        dplyr::sample_n(41) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbMinimal |>
        dplyr::filter(referenceType == "original") |>
        dplyr::sample_n(58) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
    )
}

# additional again (new process) (needs 2_validating)
if (exists("realMetaSample")) {
  openDbMinimal <- dplyr::anti_join(inhouseDbFull, realMetaSample)
}
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (exists("realMetaSample")) {
  additionalSetBis <-
    dplyr::bind_rows(
      A <- openDbMinimal |>
        dplyr::filter(referenceType == "title") |>
        dplyr::sample_n(2) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbMinimal |>
        dplyr::filter(referenceType == "publishingDetails") |>
        dplyr::sample_n(15) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbMinimal |>
        dplyr::filter(referenceType == "split") |>
        dplyr::sample_n(3) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbMinimal |>
        dplyr::filter(referenceType == "pubmed") |>
        dplyr::sample_n(10) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      E <- openDbMinimal |>
        dplyr::filter(referenceType == "doi") |>
        dplyr::sample_n(2) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        )
    ) |>
    dplyr::select(
      database,
      organismOriginal = organismValue,
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
      curator,
      validated,
      comments
    )
}

log_debug("ensuring directories exist")
create_dir(export = pathDataInterimTablesAnalyzed)

log_debug("exporting ...")
log_debug(pathDataInterimTablesAnalyzedSampleAllONPDB)
if (exists("sampleONPDB")) {
  write.table(
    x = sampleONPDB,
    file = pathDataInterimTablesAnalyzedSampleAllONPDB,
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

## 2nd iteration where some knapsack entries were needed
log_debug(pathDataInterimTablesAnalyzedSampleKnapsack)
if (exists("sampleKnapsack")) {
  write.table(
    x = sampleKnapsack,
    file = pathDataInterimTablesAnalyzedSampleKnapsack,
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

log_debug(file.path(
  pathDataInterimTablesAnalyzed,
  "samplePublishingDetails.tsv"
))
if (exists("sampleONPDB_publishingDetails")) {
  write.table(
    x = sampleONPDB_publishingDetails,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "samplePublishingDetails.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

log_debug(file.path(
  pathDataInterimTablesAnalyzed,
  "additionalSet.tsv"
))

if (exists("additionalSet")) {
  write.table(
    x = additionalSet,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "additionalSet.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("additionalSetBis")) {
  write.table(
    x = additionalSetBis,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "additionalSetBis.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

## 3rd iteration, 2_validating needs to run first.
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (exists("openDbClean2")) {
  additionalSetTer <-
    dplyr::bind_rows(
      A <- openDbClean2 |>
        dplyr::filter(referenceType == "title") |>
        dplyr::sample_n(49) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbClean2 |>
        dplyr::filter(referenceType == "publishingDetails") |>
        dplyr::sample_n(12) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbClean2 |>
        dplyr::filter(referenceType == "pubmed") |>
        dplyr::sample_n(59) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        )
    ) |>
    dplyr::select(
      database,
      organismOriginal = organismValue,
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
      curator,
      validated,
      comments
    )
}

## 4th iteration, 2_validating needs to run first.
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (exists("openDbClean3")) {
  additionalSetTetra <-
    dplyr::bind_rows(
      A <- openDbClean2 |>
        dplyr::filter(referenceType == "pubmed") |>
        dplyr::sample_n(40) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbClean2 |>
        dplyr::filter(referenceType == "split") |>
        dplyr::sample_n(220) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbClean2 |>
        dplyr::filter(referenceType == "original") |>
        dplyr::sample_n(25) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbClean2 |>
        dplyr::filter(referenceType == "title") |>
        dplyr::sample_n(25) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbClean2 |>
        dplyr::filter(referenceType == "publishingDetails") |>
        dplyr::sample_n(50) |>
        dplyr::mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        )
    ) |>
    dplyr::select(
      database,
      organismOriginal = organismValue,
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
      curator,
      validated,
      comments
    )
}

if (exists("additionalSetTer")) {
  write.table(
    x = additionalSetTer,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "additionalSetTer.tsv"
    ),
    na = "",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

if (exists("additionalSetTetra")) {
  write.table(
    x = additionalSetTetra,
    file = file.path(
      pathDataInterimTablesAnalyzed,
      "additionalSetTetra.tsv"
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
