# title: "Knapsack scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(readr)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(xml2)

source("r/progressr.R")

# get paths
database <- databases$get("knapsack")

url <-
  "http://www.knapsackfamily.com/knapsack_core/information.jsp?word=C00"

xs <- 1:63574

GetKnapSack <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          cd_id <- stringr::str_pad(
            string = x,
            width = 6,
            pad = "0"
          )
          url_id <- paste(url, cd_id) |>
            gsub(pattern = "\\s", replacement = "")
          df0 <- rvest::read_html(url_id) |>
            rvest::html_element(xpath = "/html/body/div/div[2]/div[2]/table") |>
            xml2::xml_child(2) |>
            xml2::xml_child(1) |>
            xml2::xml_child(1)

          df1 <- df0 |>
            rvest::html_table(fill = TRUE) |>
            head(8) |>
            dplyr::select(1, 2) |>
            t()

          colnames(df1) <- df1[1, ]

          df1 <- data.frame(df1) |>
            dplyr::slice(2) |>
            dplyr::mutate(joincol = url_id)

          df2 <- df0 |>
            xml2::xml_child(10) |>
            xml2::xml_child(2) |>
            xml2::xml_child(1) |>
            rvest::html_table(fill = TRUE)

          i <- 0:(nrow(df2) - 1) |>
            unlist() |>
            data.frame() |>
            dplyr::mutate(joincol = url_id) |>
            dplyr::select(i = 1, joincol)

          df_for_ref <- dplyr::full_join(df1, i) |>
            dplyr::mutate(link = paste0(joincol, "&key=", i)) |>
            dplyr::select(-joincol, -i)

          return(df_for_ref)
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}

GetKnapSackRef <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          df3 <- rvest::read_html(x) |>
            rvest::html_element(xpath = "/html/body/div/div[2]/div[2]/table") |>
            xml2::xml_child(2) |>
            xml2::xml_child(2) |>
            xml2::xml_child(7) |>
            rvest::html_table(fill = TRUE) |>
            t()

          colnames(df3) <- df3[1, ]

          df3 <- data.frame(df3) |>
            dplyr::slice(2)

          df4 <- cbind(c("link" = x, df3)) |>
            t() |>
            data.frame() |>
            dplyr::mutate_all(as.character)

          return(df4)
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}

df1 <- GetKnapSack(xs = xs) |>
  progressr::with_progress()

KnapSackTable <- dplyr::bind_rows(df1[!is.na(df1)])

xs <- KnapSackTable$link

df3 <- GetKnapSackRef(xs = xs) |>
  progressr::with_progress()

df4 <- bind_rows(df3[!is.na(df3)])

KNApSAcK_db <- KnapSackTable |>
  dplyr::full_join(df4)

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, KNApSAcK_db)
