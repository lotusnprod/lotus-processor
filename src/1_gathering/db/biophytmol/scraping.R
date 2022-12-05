# title: "BIOPHYTMOL scrapeR"

# loading paths
source("paths.R")

library(data.table)
library(dplyr)
library(rvest) # provides read_html
library(tidyr) # provides pivot_wider
library(xml2)

# get paths
database <- databases$get("biophytmol")

url <-
  "http://crdd.osdd.net/servers/biophytmol/search-biophytmol.php?compound_id="

X <- (1001:4154)

getbiophyt <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "&type=compound_id")
    url_id <- gsub(
      pattern = "\\s",
      replacement = "",
      x = url_id
    )
    df1 <- rvest::read_html(x = url_id) |>
      rvest::html_element("body") |>
      xml2::xml_child("table[3]") |>
      rvest::html_table(fill = TRUE)
  })
}

BIOPHYTMOL <- invisible(lapply(
  FUN = getbiophyt,
  X = X
))

BIOPHYTMOL_2 <- BIOPHYTMOL[BIOPHYTMOL != "Timed out!"]

BIOPHYTMOL_3 <- dplyr::bind_rows(BIOPHYTMOL_2, .id = "column_label")

BIOPHYTMOL_4 <- BIOPHYTMOL_3 |>
  dplyr::select(1:3) |>
  dplyr::filter(!is.na(X2)) |>
  dplyr::group_by(column_label) |>
  tidyr::pivot_wider(
    names_from = X1,
    values_from = X2
  ) |>
  dplyr::filter(!is.na(SMILES)) |>
  dplyr::select(
    uniqueid = `Compound ID`,
    name = `Active Compound Identified`,
    biologicalsource = `Plant Source`,
    biologicalpart = `Plant Part Used`,
    extract = Extract,
    pubchem = `PubChem ID`,
    smiles = SMILES,
    pubmed = `PubMed ID [Source Literature]`,
    reference = `Reference(s)`
  ) |>
  dplyr::mutate(reference = paste(pubmed, reference, sep = "§")) |>
  dplyr::select(
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
create_dir(export = database$sourceFiles$tsv)
database$writeFile(database$sourceFiles$tsv, BIOPHYTMOL_4)
