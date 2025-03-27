# title: "PHYTOHUB scrapeR"

# loading paths
source("paths.R")
source("r/y_as_na.R")

library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html

source("r/progressr.R")

# get paths
database <- databases$get("phytohub")

url <- "https://phytohub.eu/entries/PHUB"

xs <- 1:2527

getphytohub <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
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
          url_id <- paste(url, cd_id)
          url_id <- gsub(
            pattern = "\\s",
            replacement = "",
            x = url_id
          )
          sample <- rvest::read_html(url_id)
          scrape1 <-
            rvest::html_element(
              sample,
              xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[2]"
            ) |>
            rvest::html_text()
          scrape2 <-
            rvest::html_element(
              sample,
              xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[11]/pre"
            ) |>
            rvest::html_text()
          scrape3 <-
            rvest::html_element(
              sample,
              xpath = "//*[@id=\"fs\"]/div/div/div/table"
            ) |>
            rvest::html_table()
          scrape4 <- scrape3[, 2]
          scrape5 <-
            rvest::html_element(
              sample,
              xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[12]/pre"
            ) |>
            rvest::html_text()

          df <- cbind(scrape1, scrape2, scrape5, scrape4)
          final_df <- data.frame(df)
          return(final_df)
        },
        error = function(e) {
          "Timed out!"
        }
      )
    }
  )
}

PHYTOHUB <- getphytohub(xs = xs)

PHYTOHUB_2 <- PHYTOHUB[PHYTOHUB != "Timed out!"]

PHYTOHUB_3 <- dplyr::bind_rows(PHYTOHUB_2)

PHYTOHUB_4 <- PHYTOHUB_3 |>
  dplyr::select(
    name = scrape1,
    inchi = scrape2,
    smiles = scrape5,
    biologicalsource = Name
    # name_precursor = Precursor,
    # biologicalsource_precursor = Food.Source
  ) |>
  # dplyr::filter(is.na(name_precursor)) |>
  dplyr::select(
    name,
    inchi,
    smiles,
    biologicalsource
  )

url <- "http://phytohub.eu/entry_food_sources/"

xs <- 1:3257

getphytohubref <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
          p(sprintf("x=%g", x))
          cd_id <- x
          url_id <- paste(url, cd_id)
          url_id <- gsub(
            pattern = "\\s",
            replacement = "",
            x = url_id
          )
          sample <- rvest::read_html(url_id)
          scrape1 <-
            rvest::html_elements(sample, xpath = "/html/body/main/div/h1") %>%
            rvest::html_text()
          scrape2 <-
            rvest::html_elements(
              sample,
              xpath = "/html/body/main/blockquote"
            ) %>%
            rvest::html_text()

          df <- cbind(scrape1, scrape2)

          final_df <- as.data.frame(df)

          return(final_df)
        },
        error = function(e) {
          "Timed out!"
        }
      )
    }
  )
}

PHYTOHUB_REF <- getphytohubref(xs = xs)

PHYTOHUB_5 <- PHYTOHUB_REF[PHYTOHUB_REF != "Timed out!"]

PHYTOHUB_6 <- dplyr::bind_rows(PHYTOHUB_5)

colnames(PHYTOHUB_6) <- c("pair", "reference")

PHYTOHUB_7 <- PHYTOHUB_6 |>
  dplyr::mutate(
    joining_col = gsub(
      pattern = "Publications for ",
      replacement = "",
      x = pair
    )
  ) |>
  dplyr::mutate(
    joining_col = gsub(
      pattern = "being present in ",
      replacement = "",
      x = joining_col
    )
  )

PHYTOHUB_8 <- PHYTOHUB_4 |>
  dplyr::mutate(joining_col = paste(name, biologicalsource, sep = " "))

PHYTOHUB_9 <- dplyr::full_join(PHYTOHUB_7, PHYTOHUB_8) |>
  dplyr::select(name, inchi, smiles, biologicalsource, reference) |>
  dplyr::filter(!is.na(name) | !is.na(inchi) | !is.na(smiles))

PHYTOHUB_9$reference <- y_as_na(PHYTOHUB_9$reference, "")

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, PHYTOHUB_9)
