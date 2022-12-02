# title: "wakankensaku scrapeR"

## NOT WORKING ANYMORE

# loading paths
source("paths.R")

library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("wakankensaku")

url <-
  "https://wakankensaku.inm.u-toyama.ac.jp/wiki/Persist:CompoundList"

WAKANKENSAKU <- rvest::read_html(url) |>
  xml2::xml_child(2) |>
  xml2::xml_child(1) |>
  xml2::xml_child(1) |>
  xml2::xml_child(1) |>
  xml2::xml_child(3) |>
  xml2::xml_child(4) |>
  xml2::xml_child(1) |>
  rvest::html_table(fill = TRUE)

WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub(
      pattern = "\r\n",
      replacement = " ",
      x = x
    )
  })
WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub(
      pattern = "\r",
      replacement = " ",
      x = x
    )
  })
WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub(
      pattern = "\n",
      replacement = " ",
      x = x
    )
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, WAKANKENSAKU)
