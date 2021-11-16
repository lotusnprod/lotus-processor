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
  read_delim(
    file = pathDataInterimTablesCuratedTable,
    delim = "\t",
    col_types = cols(.default = "c")
  )

log_debug("sampling ...")
log_debug("... DOI")
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (nrow(openDbMinimal %>%
  filter(referenceType == "doi")) >= 30) {
  sampleONPDB_doi <- openDbMinimal %>%
    filter(referenceType == "doi") %>%
    sample_n(30) %>%
    mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "original")) >= 30) {
  sampleONPDB_original <- openDbMinimal %>%
    filter(referenceType == "original") %>%
    sample_n(30) %>%
    mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "pubmed")) >= 30) {
  sampleONPDB_pubmed <- openDbMinimal %>%
    filter(referenceType == "pubmed") %>%
    sample_n(30) %>%
    mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "split")) >= 30) {
  sampleONPDB_split <- openDbMinimal %>%
    filter(referenceType == "split") %>%
    sample_n(30) %>%
    mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "title")) >= 30) {
  sampleONPDB_title <- openDbMinimal %>%
    filter(referenceType == "title") %>%
    sample_n(30) %>%
    mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "publishingDetails")) >= 30) {
  sampleONPDB_publishingDetails <- openDbMinimal %>%
    filter(referenceType == "publishingDetails") %>%
    sample_n(30) %>%
    mutate(
      curator = "AR",
      validated = NA,
      comments = NA
    )
}

# get0 if to avoid error in minimal mode if df not present
sampleONPDB <- bind_rows(
  get0(x = "sampleONPDB_doi"),
  get0(x = "sampleONPDB_original"),
  get0(x = "sampleONPDB_publishingDetails"),
  get0(x = "sampleONPDB_pubmed"),
  get0(x = "sampleONPDB_split"),
  get0(x = "sampleONPDB_title")
) %>%
  mutate_all(as.character)

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
#   filter(database == "kna_1") %>%
#   sample_n(150) %>%
#   mutate(
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
if (nrow(openDbMinimal %>%
  filter(referenceType == "title")) >= 85 &
  nrow(openDbMinimal %>%
    filter(referenceType == "publishingDetails")) >= 25 &
  nrow(openDbMinimal %>%
    filter(referenceType == "split")) >= 41 &
  nrow(openDbMinimal %>%
    filter(referenceType == "publishingDetails")) >= 58) {
  additionalSet <-
    bind_rows(
      A <- openDbMinimal %>%
        filter(referenceType == "title") %>%
        sample_n(85) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbMinimal %>%
        filter(referenceType == "publishingDetails") %>%
        sample_n(25) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbMinimal %>%
        filter(referenceType == "split") %>%
        sample_n(41) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbMinimal %>%
        filter(referenceType == "original") %>%
        sample_n(58) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
    )
}

# additional again (new process) (needs 2_validating)
if (exists("realMetaSample")) {
  openDbMinimal <- anti_join(inhouseDbFull, realMetaSample)
}
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (exists("realMetaSample")) {
  additionalSetBis <-
    bind_rows(
      A <- openDbMinimal %>%
        filter(referenceType == "title") %>%
        sample_n(2) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbMinimal %>%
        filter(referenceType == "publishingDetails") %>%
        sample_n(15) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbMinimal %>%
        filter(referenceType == "split") %>%
        sample_n(3) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      D <- openDbMinimal %>%
        filter(referenceType == "pubmed") %>%
        sample_n(10) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      E <- openDbMinimal %>%
        filter(referenceType == "doi") %>%
        sample_n(2) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
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
ifelse(
  test = !dir.exists(pathDataInterimTablesAnalyzed),
  yes = dir.create(pathDataInterimTablesAnalyzed),
  no = paste(pathDataInterimTablesAnalyzed, "exists")
)

log_debug("exporting ...")
log_debug(pathDataInterimTablesAnalyzedSampleAllONPDB)
if (exists("sampleONPDB")) {
  write.table(
    x = sampleONPDB,
    file = pathDataInterimTablesAnalyzedSampleAllONPDB,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

log_debug(pathDataInterimTablesAnalyzedSampleKnapsack)
if (exists("sampleKnapsack")) {
  write.table(
    x = sampleKnapsack,
    file = pathDataInterimTablesAnalyzedSampleKnapsack,
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
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}

# additional again
set.seed(
  seed = 42,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion"
)
if (exists("openDbClean2")) {
  additionalSetTer <-
    bind_rows(
      A <- openDbClean2 %>%
        filter(referenceType == "title") %>%
        sample_n(49) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      B <- openDbClean2 %>%
        filter(referenceType == "publishingDetails") %>%
        sample_n(12) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
        ),
      C <- openDbClean2 %>%
        filter(referenceType == "pubmed") %>%
        sample_n(59) %>%
        mutate(
          curator = "AR",
          validated = NA,
          comments = NA
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
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}


end <- Sys.time()

log_debug("Script finished in", format(end - start))
