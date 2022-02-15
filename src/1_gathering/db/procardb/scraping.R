# title: "PROCARDB scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(tidyr) # provides pivot_wider

# get paths
database <- databases$get("procardb")

url <-
  "http://bioinfo.imtech.res.in/servers/procardb/?c=carotenoides&m=getDetail&id=C"

url_bio <-
  "http://bioinfo.imtech.res.in/servers/procardb/?c=organisms&m=getDetail&id=R"

X <- (1:1942)
X_bio <- (1:403)

getprocardb <- function(X) {
  tryCatch(
    {
      cd_id <- str_pad(X, 4, pad = "0")
      url_id <- paste(url, cd_id)
      url_id <- gsub("\\s", "", url_id)
      df1 <- read_html(url_id) %>%
        html_element(xpath = "/html/body/div[3]/div/div/table") %>%
        html_table(., fill = TRUE)

      return(df1)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

PROCARDB <- invisible(
  pbmclapply(
    FUN = getprocardb,
    X = X,
    mc.cores = numCores,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

PROCARDB_2 <- PROCARDB[PROCARDB != "Timed out!"]

for (i in seq_along(PROCARDB_2)) {
  colnames(PROCARDB_2[[i]]) <-
    c("CAROTENOID INFO", "CAROTENOID INFO 2")
}

PROCARDB_2 <- bind_rows(PROCARDB_2, .id = "column_label")

PROCARDB_2 <- PROCARDB_2 %>%
  group_by(column_label) %>%
  pivot_wider(
    names_from = 2,
    values_from = 3
  )

urls_1 <- X
ids_1 <- PROCARDB_2$column_label
list_1 <- urls_1[which(!urls_1 %in% ids_1)]

getprocardb_bio <- function(X_bio) {
  tryCatch(
    {
      cd_id_bio <- str_pad(X_bio, 3, pad = "0")
      url_id_bio <- paste(url_bio, cd_id_bio)
      url_id_bio <- gsub("\\s", "", url_id_bio)
      df1 <- read_html(url_id_bio) %>%
        html_element(xpath = "/html/body/div[3]/div/div/div[2]/table") %>%
        html_table(., fill = TRUE)

      return(df1)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

PROCARDB_3 <- invisible(
  pbmclapply(
    FUN = getprocardb_bio,
    X = X_bio,
    mc.cores = numCores,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

PROCARDB_4 <- PROCARDB_3[PROCARDB_3 != "Timed out!"]

PROCARDB_4 <- bind_rows(PROCARDB_4, .id = "column_label")

getprocardb_names <- function(X_bio) {
  tryCatch(
    {
      cd_id_bio <- str_pad(X_bio, 3, pad = "0")
      url_id_bio <- paste(url_bio, cd_id_bio)
      url_id_bio <- gsub("\\s", "", url_id_bio)
      df1 <- read_html(url_id_bio) %>%
        html_element(xpath = "/html/body/div[3]/div/div/div") %>%
        html_text()

      return(df1)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

PROCARDB_5 <- invisible(
  pbmclapply(
    FUN = getprocardb_names,
    X = X_bio,
    mc.cores = numCores,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

PROCARDB_6 <- PROCARDB_5[!is.na(PROCARDB_5)]

PROCARDB_6 <-
  data.frame(matrix(
    unlist(PROCARDB_6),
    nrow = length(PROCARDB_6),
    byrow = TRUE
  ))

PROCARDB_6 <- PROCARDB_6 %>%
  select(biologicalsource = 1) %>%
  mutate(column_label = row.names(.))

PROCARDB_final <- full_join(PROCARDB_6, PROCARDB_4) %>%
  select(
    biologicalsource,
    `CAROTENOID NAME`
  )

PROCARDB_final <- full_join(PROCARDB_final, PROCARDB_2) %>%
  distinct(biologicalsource, `CAROTENOID NAME`, .keep_all = TRUE)

PROCARDB_final[] <-
  lapply(PROCARDB_final, function(x) {
    gsub("\r\n", " ", x)
  })
PROCARDB_final[] <-
  lapply(PROCARDB_final, function(x) {
    gsub("\r", " ", x)
  })
PROCARDB_final[] <-
  lapply(PROCARDB_final, function(x) {
    gsub("\n", " ", x)
  })
PROCARDB_final[] <-
  lapply(PROCARDB_final, function(x) {
    gsub("\t", " ", x)
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, PROCARDB_final)
