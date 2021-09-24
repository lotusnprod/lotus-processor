# title: "CarotenoidDB scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")
source("r/y_as_na.R")

library(data.table)
library(dplyr)
library(parallel)
library(pbmcapply)
library(rvest) # provides read_html
library(tidyr) # provides pivot_wider
library(xml2)

# get paths
database <- databases$get("carotenoiddb")

## files
ids <- read_delim(
  file = database$sourceFiles$tsvInchi,
  col_names = FALSE,
) %>%
  mutate_all(as.character)

url <- "http://carotenoiddb.jp/Entries/"

X <- ids$X1

getcarotenoid <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, ".html")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_element("body") %>%
      html_element("div") %>%
      html_element("td.fr2") %>%
      xml_child(1) %>%
      html_table(fill = TRUE)
  })
}

CAROTENOIDDB <- invisible(
  pbmclapply(
    FUN = getcarotenoid,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 1),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

CAROTENOIDDB_2 <- bind_rows(CAROTENOIDDB) %>%
  select(X1, X2)

CAROTENOIDDB_2$level <-
  as.numeric(gl(nrow(CAROTENOIDDB_2) / 31, 31))

colnames(CAROTENOIDDB_2) <- c("name", "value", "level")

CAROTENOIDDB_2$name <- y_as_na(CAROTENOIDDB_2$name, "")
CAROTENOIDDB_2$value <- y_as_na(CAROTENOIDDB_2$value, "")

CAROTENOIDDB_3 <- CAROTENOIDDB_2 %>%
  filter(!grepl("^CA0", name)) %>%
  group_by(level) %>%
  pivot_wider(
    names_from = name,
    values_from = value
  ) %>%
  ungroup() %>%
  select(-level)

CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub("\r\n", " ", x)
  })
CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub("\r", " ", x)
  })
CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub("\n", " ", x)
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, CAROTENOIDDB_3)