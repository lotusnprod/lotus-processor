cat("This script compares OTL (Open Tree of Life) IDs obtained via gnverify and via rotl API \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... functions \n")
source("r/log.R")
source("r/vroom_safe.R")

cat("loading ... \n")
cat("... libraries \n")
library(DBI)
library(rotl)
library(RSQLite)
library(tidyverse)

canonical_name_colname <- "organismCleaned"

fullOrganisms <-
  vroom_read_safe(path = pathDataInterimDictionariesOrganismDictionary)

drv <- SQLite()

db <- dbConnect(
  drv = drv,
  dbname = pathDataInterimDictionariesOrganismDictionaryOTL
)

otlVersion <- dbGetQuery(
  conn = db,
  statement = "SELECT
  taxa_names.canonical_name AS organismCleaned,
  taxa_otl.ott_id AS organismCleanedId
  FROM taxa_names
  LEFT JOIN taxa_otl
  ON taxa_names.search_string = taxa_otl.search_string"
) %>%
  filter(organismCleaned %in% fullOrganisms$organismCleaned) %>%
  filter(!is.na(organismCleanedId)) %>%
  distinct()

cat(nrow(otlVersion), "for rotl API version \n")

gnverifyVersion <- fullOrganisms %>%
  filter(organismCleaned_dbTaxo == "Open Tree of Life") %>%
  filter(!is.na(organismCleaned_id)) %>%
  distinct(organismCleaned, organismCleaned_id) %>%
  mutate(organismCleanedId = as.integer(organismCleaned_id))

cat(nrow(gnverifyVersion), "for gnverify version \n")

diff <- anti_join(otlVersion, gnverifyVersion)

cat("so it seems that rotl version has \n", nrow(diff), "more results \n")

cat("but then if we have a closer look and only go for distinct IDs... \n")
cat("just an example to show what is meant: \n")
cat("applying filter(organismCleanedId == \"65272\") to both tables \n")

otlStrepto <- otlVersion %>%
  filter(organismCleanedId == "65272")

cat("rotl table \n")
otlStrepto

gnverifyStrepto <- gnverifyVersion %>%
  filter(organismCleanedId == "65272")

cat("gnverify table \n")
gnverifyStrepto

cat("so when keeping only one name per ID... \n")

otlVersionDistinct <- otlVersion %>%
  left_join(., gnverifyVersion %>% mutate(isInBoth = "Y")) %>% ## to keep the same synonym for both
  arrange(isInBoth) %>%
  distinct(organismCleanedId, .keep_all = TRUE)
cat(nrow(otlVersionDistinct), "distinct IDs for rotl API version \n")

gnverifyVersionDistinct <- gnverifyVersion %>%
  distinct(organismCleanedId, .keep_all = TRUE)

cat(
  nrow(gnverifyVersionDistinct),
  "distinct IDs for gnverify version \n"
)

diffDistinct_1 <-
  anti_join(otlVersionDistinct, gnverifyVersionDistinct)

diffDistinct_2 <-
  anti_join(gnverifyVersionDistinct, otlVersionDistinct)

cat(
  "so it seems that rotl version has \n",
  nrow(diffDistinct_1),
  "names that are not in gnverify \n"
)
cat(
  "and that gnverify version has \n",
  nrow(diffDistinct_2),
  "names that are not in rotl \n"
)

cat("Conclusion: almost the same but very interesting for synonyms filtering \n")

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
