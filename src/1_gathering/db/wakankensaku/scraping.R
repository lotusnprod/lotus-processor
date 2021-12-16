# title: "wakankensaku scrapeR"

## NOT WORKING ANYMORE

# loading paths
source("paths.R")
source("r/parallel.R")

library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("wakankensaku")

url <-
  "https://wakankensaku.inm.u-toyama.ac.jp/wiki/Persist:CompoundList"

WAKANKENSAKU <- read_html(url) %>%
  xml_child(2) %>%
  xml_child(1) %>%
  xml_child(1) %>%
  xml_child(1) %>%
  xml_child(3) %>%
  xml_child(4) %>%
  xml_child(1) %>%
  html_table(fill = TRUE)

WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub("\r\n", " ", x)
  })
WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub("\r", " ", x)
  })
WAKANKENSAKU[] <-
  lapply(WAKANKENSAKU, function(x) {
    gsub("\n", " ", x)
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, WAKANKENSAKU)
