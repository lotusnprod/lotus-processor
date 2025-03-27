# title: "SWMD scrapeR"

# loading paths
source("paths.R")
source("r/y_as_na.R")

library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

source("r/progressr.R")

# get paths
database <- databases$get("swmd")

url <- "http://www.swmd.co.in/search.php?No="

xs <- list.files(path = pathDataExternalDbSourceSwmdDirectory) |>
  gsub(
    pattern = ".mol",
    replacement = ""
  )

getswmd <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", as.numeric(x))) ## little hack
          cd_id <- x
          url_id <- paste0(url, cd_id)

          df0 <- rvest::read_html(url_id) |>
            rvest::html_element(xpath = "body")

          df1 <- df0 |>
            rvest::html_element(xpath = "div[1]/table[3]") |>
            rvest::html_table(fill = TRUE) |>
            dplyr::select(1, 2)

          df2 <- df0 |>
            rvest::html_element(xpath = "div[2]/table") |>
            rvest::html_table(fill = TRUE)

          df3 <- df0 |>
            rvest::html_element(xpath = "div[4]/div/table") |>
            rvest::html_table(fill = TRUE)

          df4 <- rbind(df1, df2, df3)

          return(df4)
        },
        error = function(e) {
          "Timed out!"
        }
      )
    }
  )
}

SWMD <- getswmd(xs = xs)

SWMD <- SWMD[SWMD != "Timed out!"]

SWMD_2 <- dplyr::bind_rows(SWMD)

SWMD_2$level <- as.numeric(gl(nrow(SWMD_2) / 21, 21))

colnames(SWMD_2) <- c("name", "value", "level")

SWMD_2$name <- y_as_na(SWMD_2$name, "")
SWMD_2$value <- y_as_na(SWMD_2$value, "")

SWMD_3 <- SWMD_2 |>
  dplyr::filter(
    !stringr::str_detect(
      string = name,
      pattern = "Accession Number\r\n"
    )
  ) |>
  dplyr::filter(!is.na(name)) |>
  dplyr::group_by(level) |>
  tidyr::pivot_wider(
    names_from = name,
    values_from = value
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-level)

SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub(
    pattern = "\r\n",
    replacement = " ",
    x = x
  )
})
SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub(
    pattern = "\r",
    replacement = " ",
    x = x
  )
})
SWMD_3[] <- lapply(SWMD_3, function(x) {
  gsub(
    pattern = "\n",
    replacement = " ",
    x = x
  )
})

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, SWMD_3)
