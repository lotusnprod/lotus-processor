# title: "SANCDB scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

# get paths
database <- databases$get("sancdb")

url <- "https://sancdb.rubi.ru.ac.za/compounds/"

X <- (1:1000)

getsanc <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "/")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_element("body") %>%
      html_element("div#wrap") %>%
      html_element("div#content.content") %>%
      html_element("div#pt-main.pt-perspective") %>%
      html_text()
  })
}

SANCDB <- invisible(
  lapply(
    FUN = getsanc,
    X = X
  )
) %>%
  as.data.table() %>%
  t() %>%
  cSplit("V1", "\n")

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, SANCDB)
