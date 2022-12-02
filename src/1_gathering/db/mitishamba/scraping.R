# title: "MITISHAMBA scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("mitishamba")

X <- 1:1102

url <- "http://mitishamba.uonbi.ac.ke/details.php?mol=00"

GetMitishamba <- function(X) {
  tryCatch(
    {
      cd_id <- str_pad(
        string = X,
        width = 6,
        pad = "0"
      )
      url_id <- paste(url, cd_id)
      url_id <- gsub(pattern = "\\s", replacement = "", url_id)
      df1 <- rvest::read_html(url_id) |>
        rvest::html_element("body") |>
        rvest::html_element(xpath = '//*[@id="wrapper-search"]') |>
        xml2::xml_child(5) |>
        rvest::html_table()

      df2 <- t(df1)
      colnames(df2) <- df2[1, ]
      df3 <- data.frame(df2) %>%
        filter(rownames(.) == "X2")
      df3[setdiff(
        row(df3),
        c(
          "mw",
          "mmff",
          "logp",
          "psa",
          "smiles",
          "rotatable_bonds",
          "hydrogen_acceptors",
          "hydrogen_donors",
          "heavy_atoms",
          "plant_family",
          "plant_species",
          "plant_part",
          "compound_type",
          "common_name",
          "authors",
          "url",
          "name",
          "place_of_collection"
        )
      )] <- NA
      return(df3)
    },
    error = function(e) {
      "Timed out!"
    }
  )
}

MITISHAMBA <- invisible(lapply(
  FUN = GetMitishamba,
  X = X
))

MITISHAMBA_2 <- MITISHAMBA[MITISHAMBA != "Timed out!"]

MITISHAMBA_3 <- dplyr::bind_rows(MITISHAMBA_2)

MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub(pattern = "\r\n", replacement = " ", x = x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub(pattern = "\r", replacement = " ", x = x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub(pattern = "\n", replacement = " ", x = x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub(pattern = "\t", replacement = " ", x = x)
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, MITISHAMBA_3)
