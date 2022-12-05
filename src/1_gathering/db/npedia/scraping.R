# title: "NPEDIA scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(tidyr) # provides pivot_wider
library(xml2)

# get paths
database <- databases$get("npedia")

url <- "http://www.cbrg.riken.jp/npedia/details.php?ID="

X <- 1:83797

getnpedia <- function(X) {
  tryCatch(
    {
      cd_id <- stringr::str_pad(
        string = X,
        width = 5,
        pad = "0"
      )
      url_id <- paste(url, cd_id, "&TAB=Basic")
      url_id <- gsub(pattern = "\\s", replacement = "", url_id)
      df1 <- rvest::read_html(url_id) |>
        rvest::html_element(xpath = "/html/body/table[2]") |>
        rvest::html_table(fill = TRUE)

      return(df1)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

NPEDIA <- invisible(lapply(
  FUN = getnpedia,
  X = X
))

NPEDIA_2 <- NPEDIA[NPEDIA != "Timed out!"]

NPEDIA_2 <- dplyr::bind_rows(NPEDIA_2, .id = "column_label")

NPEDIA_2$X1 <-
  gsub(
    pattern = "=.*",
    replacement = "",
    x = NPEDIA_2$X1
  )

NPEDIA_2 <- NPEDIA_2 |>
  dplyr::group_by(column_label) |>
  tidyr::pivot_wider(
    names_from = X1,
    values_from = X2
  ) |>
  dplyr::ungroup()

urls_1 <- stringr::str_pad(
  string = X,
  width = 5,
  pad = "0"
)
ids_1 <- NPEDIA_2$ID
list_1 <- urls_1[which(!urls_1 %in% NPEDIA_2$ID)]

getnpedia_2 <- function(X) {
  tryCatch(
    {
      cd_id <- stringr::str_pad(
        string = X,
        width = 5,
        pad = "0"
      )
      url_id <- paste(url, cd_id, "&TAB=Origin")
      url_id <- gsub(pattern = "\\s", replacement = "", url_id)
      df1 <- rvest::read_html(url_id) |>
        rvest::html_element(xpath = "/html/body/table[2]") |>
        rvest::html_table(fill = TRUE)

      return(df1)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

NPEDIA_3 <- invisible(lapply(
  FUN = getnpedia_2,
  X = X
))

NPEDIA_4 <- NPEDIA_3[NPEDIA_3 != "Timed out!"]

NPEDIA_4 <- dplyr::bind_rows(NPEDIA_4, .id = "column_label")

NPEDIA_4 <- NPEDIA_4 |>
  dplyr::group_by(column_label) |>
  tidyr::pivot_wider(
    names_from = X1,
    values_from = X2
  ) |>
  dplyr::ungroup() |>
  tidyr::unnest()

urls_2 <- X
ids_2 <- NPEDIA_4$column_label
list_2 <- urls_2[which(!urls_2 %in% NPEDIA_4$column_label)]

NPEDIA_4 <- NPEDIA_4 |>
  dplyr::filter(!is.na(Source)) |>
  dplyr::filter(Source != "")

NPEDIA_final <- dplyr::full_join(NPEDIA_2, NPEDIA_4)

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, NPEDIA_final)
