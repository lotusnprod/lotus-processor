# title: "TMDB scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

# get paths
database <- databases$get("tmdb")

url <- "http://pcsb.ahau.edu.cn:8080/TCDB/f/browseDetail?id="

X <- 1:1473

gettmdb <- function(X) {
  tryCatch(
    {
      cd_id <- X
      url_id <- paste(url, cd_id)
      url_id <- gsub(
        pattern = "\\s",
        replacement = "",
        x = url_id
      )
      sample <- read_html(url_id)
      scrape1 <-
        rvest::html_elements(sample, xpath = "/html/body/div[1]/div/table") |>
        rvest::html_table(fill = TRUE)

      scrape2 <- scrape1[[1]]
      return(scrape2)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

TMDB <- invisible(lapply(
  FUN = gettmdb,
  X = X
))

TMDB_2 <- TMDB[TMDB != "Timed out!"]

TMDB_3 <- dplyr::bind_rows(TMDB_2)

TMDB_4 <- TMDB_3 |>
  dplyr::filter(!is.na(X1))

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, TMDB_4)
