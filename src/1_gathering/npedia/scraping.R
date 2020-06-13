# title: "NPEDIA scrapeR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

url <- 'http://www.cbrg.riken.jp/npedia/details.php?ID='

X <- (1:83797)

getnpedia <- function(X)
{
  tryCatch({
    cd_id <-  str_pad (X, 5, pad = "0")
    url_id <- paste(url, cd_id, "&TAB=Basic")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node(xpath = "/html/body/table[2]") %>%
      html_table(., fill = TRUE)
  },
  error = function(e) {
    "Timed out!"
  })
}

NPEDIA <- invisible(
  pbmclapply(
    FUN = getnpedia,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

NPEDIA_2 <- NPEDIA[NPEDIA != "Timed out!"]

NPEDIA_2 <- bind_rows(NPEDIA_2, .id = "column_label")

NPEDIA_2$X1 <- gsub("=.*", "", NPEDIA_2$X1)

NPEDIA_2 <- NPEDIA_2 %>%
  group_by(column_label) %>%
  pivot_wider(names_from = X1,
              values_from = X2) %>%
  ungroup()

urls_1 <- str_pad(X, 5, pad = "0")
ids_1 <- NPEDIA_2$ID
list_1 <- urls_1[which(!urls_1 %in% NPEDIA_2$ID)]

getnpedia_2 <- function(X)
{
  tryCatch({
    cd_id <-  str_pad (X, 5, pad = "0")
    url_id <- paste(url, cd_id, "&TAB=Origin")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node(xpath = "/html/body/table[2]") %>%
      html_table(., fill = TRUE)
  },
  error = function(e) {
    "Timed out!"
  })
}

NPEDIA_3 <- invisible(
  pbmclapply(
    FUN = getnpedia_2,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

NPEDIA_4 <- NPEDIA_3[NPEDIA_3 != "Timed out!"]

NPEDIA_4 <- bind_rows(NPEDIA_4, .id = "column_label")

NPEDIA_4 <- NPEDIA_4 %>%
  group_by(column_label) %>%
  pivot_wider(names_from = X1,
              values_from = X2) %>%
  ungroup() %>%
  unnest()

urls_2 <- X
ids_2 <- NPEDIA_4$column_label
list_2 <- urls_2[which(!urls_2 %in% NPEDIA_4$column_label)]

NPEDIA_4 <- NPEDIA_4 %>%
  filter(!is.na(Source)) %>%
  filter(Source != "")

NPEDIA_final <- full_join(NPEDIA_2, NPEDIA_4)

# exporting
write.table(
  x = NPEDIA_final,
  file = gzfile(
    description = pathDataExternalDbSourceNpediaOriginal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
