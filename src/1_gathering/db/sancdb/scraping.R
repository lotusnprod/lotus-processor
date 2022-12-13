# title: "SANCDB scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

source("r/progressr.R")

# get paths
database <- databases$get("sancdb")

url <- "https://sancdb.rubi.ru.ac.za/compounds/"

xs <- 1:1134

getsanc <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch({
        p(sprintf("x=%g", x))
        cd_id <- x
        url_id <- paste(url, cd_id, "/")
        url_id <- gsub(
          pattern = "\\s",
          replacement = "",
          x = url_id
        )
        df1 <- rvest::read_html(url_id) |>
          rvest::html_element("body") |>
          rvest::html_element("div#wrap") |>
          rvest::html_element("div#content.content") |>
          rvest::html_element("div#pt-main.pt-perspective") |>
          rvest::html_text()
      })
    }
  )
}

SANCDB <- getsanc(xs = xs) |>
  progressr::with_progress(enable = TRUE) |>
  as.data.table() |>
  t() |>
  splitstackshape::cSplit("V1", "\n")

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, SANCDB)
