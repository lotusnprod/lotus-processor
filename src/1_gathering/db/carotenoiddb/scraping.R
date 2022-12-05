# title: "CarotenoidDB scrapeR"

# loading paths
source("paths.R")
source("r/y_as_na.R")

library(data.table)
library(dplyr)
library(future)
library(future.apply)
library(progressr)
library(rvest) # provides read_html
library(tidyr) # provides pivot_wider
library(xml2)

source("r/progressr.R")

# get paths
database <- databases$get("carotenoiddb")

## files
ids <- readr::read_delim(
  file = database$sourceFiles$tsvInchi,
  col_names = FALSE
) |>
  dplyr::mutate_all(as.character)

url <- "http://carotenoiddb.jp/Entries/"

xs <- ids$X1

getcarotenoid <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch({
        p(sprintf("x=%g", as.numeric(x))) ## little hack
        cd_id <- x
        url_id <- paste(url, cd_id, ".html")
        url_id <- gsub(
          pattern = "\\s",
          replacement = "",
          x = url_id
        )
        df1 <- rvest::read_html(x = url_id) |>
          rvest::html_element("body") |>
          rvest::html_element("div") |>
          rvest::html_element("td.fr2") |>
          xml2::xml_child(1) |>
          rvest::html_table(fill = TRUE)
      })
    }
  )
}

CAROTENOIDDB <- getcarotenoid(xs = xs) |>
  progressr::with_progress()

CAROTENOIDDB_2 <- dplyr::bind_rows(CAROTENOIDDB) |>
  dplyr::select(X1, X2)

CAROTENOIDDB_2$level <-
  as.numeric(gl(nrow(CAROTENOIDDB_2) / 31, 31))

colnames(CAROTENOIDDB_2) <- c("name", "value", "level")

CAROTENOIDDB_2$name <- y_as_na(CAROTENOIDDB_2$name, "")
CAROTENOIDDB_2$value <- y_as_na(CAROTENOIDDB_2$value, "")

CAROTENOIDDB_3 <- CAROTENOIDDB_2 |>
  dplyr::filter(!grepl(pattern = "^CA0", x = name)) |>
  dplyr::group_by(level) |>
  tidyr::pivot_wider(
    names_from = name,
    values_from = value
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-level)

CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub(pattern = "\r\n", replacement = " ", x = x)
  })
CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub(pattern = "\r", replacement = " ", x = x)
  })
CAROTENOIDDB_3[] <-
  lapply(CAROTENOIDDB_3, function(x) {
    gsub(pattern = "\n", replacement = " ", x = x)
  })

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, CAROTENOIDDB_3)
