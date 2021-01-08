# title: "wakankensaku scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

groundhog.library(pbmcapply, date = groundhog.day)
library(parallel)
groundhog.library(data.table, date = groundhog.day)
groundhog.library(splitstackshape, date = groundhog.day) # provides cSplit
groundhog.library(rvest, date = groundhog.day) # provides read_html
groundhog.library(tidyverse, date = groundhog.day) # provides pivot_wider

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
database$writeFile(database$sourceFiles$tsv, WAKANKENSAKU)