# title: "ANPDB scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(xml2)

source("r/progressr.R")

# get paths
database <- databases$get("anpdb")

url <- "http://african-compounds.org/anpdb/get_compound_card/"

xs <- 1:8410

getanp <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      tryCatch(
        {
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
            xml2::xml_child(3) |>
            xml2::xml_child(1) |>
            xml2::xml_child(2) |>
            xml2::xml_child(3) |>
            xml2::xml_child(1) |>
            rvest::html_table()

          df2 <- t(df1)

          colnames(df2) <- df2[1, ]

          df3 <- data.frame(df2) %>%
            filter(rownames(.) == "x2")

          df3[setdiff(
            row(df3),
            c(
              "Image.",
              "SMILES.",
              "PubChem.",
              "Properties",
              "Source.Species.Information",
              "Predicted.toxicity.from.pkCSM",
              "Reference.information",
              "Authors.information"
            )
          )] <- NA

          return(df3)
        },
        error = function(e) {
          NA
        }
      )
    }
  )
}

ANPDB <- getanp(xs = xs) |>
  progressr::with_progress()

ANPDB_2 <- bind_rows(ANPDB[!is.na(ANPDB)])

ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "\r\n",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "\r",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "\n",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "\t",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
ANPDB_2[] <- lapply(ANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})

# exporting
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, ANPDB_2)
