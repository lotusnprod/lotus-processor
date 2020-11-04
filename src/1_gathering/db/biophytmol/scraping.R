# title: "BIOPHYTMOL scrapeR"

# loading paths
source("paths.R")
source("functions/parallel.R")

library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(tidyverse) # provides pivot_wider

# get paths
database <- databases$get("biophytmol")

url <-
  "http://crdd.osdd.net/servers/biophytmol/search-biophytmol.php?compound_id="

X <- (1001:4154)

getbiophyt <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "&type=compound_id")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node("body") %>%
      xml_child("table[3]") %>%
      html_table(., fill = TRUE)
  })
}

BIOPHYTMOL <- invisible(
  pbmclapply(
    FUN = getbiophyt,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

BIOPHYTMOL_2 <- BIOPHYTMOL[BIOPHYTMOL != "Timed out!"]

BIOPHYTMOL_3 <- bind_rows(BIOPHYTMOL_2, .id = "column_label")

BIOPHYTMOL_4 <- BIOPHYTMOL_3 %>%
  select(1:3) %>%
  filter(!is.na(X2)) %>%
  group_by(column_label) %>%
  pivot_wider(
    names_from = X1,
    values_from = X2
  ) %>%
  filter(!is.na(SMILES)) %>%
  select(
    uniqueid = `Compound ID`,
    name = `Active Compound Identified`,
    biologicalsource = `Plant Source`,
    biologicalpart = `Plant Part Used`,
    extract = Extract,
    pubchem = `PubChem ID`,
    smiles = SMILES,
    pubmed = `PubMed ID [Source Literature]`,
    reference = `Reference(s)`
  ) %>%
  mutate(reference = paste(pubmed, reference, sep = "§")) %>%
  select(
    uniqueid,
    name,
    biologicalsource,
    biologicalpart,
    extract,
    pubchem,
    smiles,
    reference
  )

# exporting
database$writeFile(database$sourceFiles$tsv, BIOPHYTMOL_4)