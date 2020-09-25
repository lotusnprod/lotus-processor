cat("This script samples some entries to then manually check their validity. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("functions.R")

cat("loading db, if running fullmode, this may take a while \n")
openDbMinimal <- read_delim(
  file = gzfile(pathDataInterimTablesCuratedTable),
  col_types = cols(.default = "c"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  data.frame()

cat("sampling ... \n")
cat("... DOI \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "doi")) >= 30)
  sampleONPDB_doi <- openDbMinimal %>%
  filter(referenceType == "doi") %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

cat("... original \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "original")) >= 30)
  sampleONPDB_original <- openDbMinimal %>%
  filter(referenceType == "original")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

cat("... PMID \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "pubmed")) >= 30)
  sampleONPDB_pubmed <- openDbMinimal %>%
  filter(referenceType == "pubmed") %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

cat("... split \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "split")) >= 30)
  sampleONPDB_split <- openDbMinimal %>%
  filter(referenceType == "split")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

cat("... title \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "title")) >= 30)
  sampleONPDB_title <- openDbMinimal %>%
  filter(referenceType == "title")  %>%
  sample_n(30) %>%
  mutate(curator = NA,
         validated = NA,
         comments = NA)

cat("... publishing details \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "publishingDetails")) >= 30)
  sampleONPDB_publishingDetails <- openDbMinimal %>%
  filter(referenceType == "publishingDetails")  %>%
  sample_n(30) %>%
  mutate(curator = "AR",
         validated = NA,
         comments = NA)

sampleONPDB <- bind_rows(
  sampleONPDB_doi,
  sampleONPDB_original,
  # sampleONPDB_publishingDetails,
  sampleONPDB_pubmed,
  sampleONPDB_split,
  sampleONPDB_title
) %>%
  mutate_all(as.character)

cat("... attributing curator \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleONPDB <- sampleONPDB[sample(nrow(sampleONPDB)),]

sampleONPDB[1:50, "curator"] <- "AR"

sampleONPDB[51:100, "curator"] <- "JB"

sampleONPDB[101:150, "curator"] <- "PMA"

cat("... knapsack entries \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
sampleKnapsack <- openDbMinimal %>%
  filter(database == "kna_1") %>%
  sample_n(150) %>%
  mutate(
    curator = sample(c("AR", "JB", "PMA"),
                     size = nrow(.),
                     replace = TRUE),
    validated = NA,
    comments = NA
  )

cat("... additional entries \n")
set.seed(seed = 42,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")
if (nrow(openDbMinimal %>%
         filter(referenceType == "title")) >= 85 &
    nrow(openDbMinimal %>%
         filter(referenceType == "publishingDetails")) >= 25 &
    nrow(openDbMinimal %>%
         filter(referenceType == "split")) >= 41 &
    nrow(openDbMinimal %>%
         filter(referenceType == "publishingDetails")) >= 58)
  additionalSet <-
  bind_rows(
    A <- openDbMinimal %>%
      filter(referenceType == "title")  %>%
      sample_n(85) %>%
      mutate(
        curator = "AR",
        validated = NA,
        comments = NA
      ),
    B <- openDbMinimal %>%
      filter(referenceType == "publishingDetails")  %>%
      sample_n(25) %>%
      mutate(
        curator = "AR",
        validated = NA,
        comments = NA
      ),
    C <- openDbMinimal %>%
      filter(referenceType == "split")  %>%
      sample_n(41) %>%
      mutate(
        curator = "AR",
        validated = NA,
        comments = NA
      ),
    D <- openDbMinimal %>%
      filter(referenceType == "original")  %>%
      sample_n(58) %>%
      mutate(
        curator = "AR",
        validated = NA,
        comments = NA
      ),
  )

cat("ensuring directories exist \n")
ifelse(
  test = !dir.exists(pathDataInterimTablesAnalysed),
  yes = dir.create(pathDataInterimTablesAnalysed),
  no = paste(pathDataInterimTablesAnalysed, "exists")
)

cat("exporting ... \n")
cat(pathDataInterimTablesAnalysedSampleAllONPDB, "\n")
write.table(
  x = sampleONPDB,
  file = pathDataInterimTablesAnalysedSampleAllONPDB,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(pathDataInterimTablesAnalysedSampleKnapsack, "\n")
write.table(
  x = sampleKnapsack,
  file = pathDataInterimTablesAnalysedSampleKnapsack,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

cat(file.path(
  pathDataInterimTablesAnalysed,
  "samplePublishingDetails.tsv"
),
"\n")
if (exists("sampleONPDB_publishingDetails"))
  write.table(
    x = sampleONPDB_publishingDetails,
    file = file.path(
      pathDataInterimTablesAnalysed,
      "samplePublishingDetails.tsv"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )

cat(file.path(pathDataInterimTablesAnalysed,
              "additionalSet.tsv"),
    "\n")
if (exists("additionalSet"))
  write.table(
    x = additionalSet,
    file = file.path(pathDataInterimTablesAnalysed,
                     "additionalSet.tsv"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
