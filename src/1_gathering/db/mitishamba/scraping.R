# title: "MITISHAMBA scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("mitishamba")

X <- (1:1102)

url <- "http://mitishamba.uonbi.ac.ke/details.php?mol=00"

GetMitishamba <- function(X) {
  tryCatch(
    {
      cd_id <- str_pad(X, 6, pad = "0")
      url_id <- paste(url, cd_id)
      url_id <- gsub("\\s", "", url_id)
      df1 <- read_html(url_id) %>%
        html_element("body") %>%
        html_element(xpath = '//*[@id="wrapper-search"]') %>%
        xml_child(5) %>%
        html_table()

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

MITISHAMBA <- invisible(
  pbmclapply(
    FUN = GetMitishamba,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 1),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

MITISHAMBA_2 <- MITISHAMBA[MITISHAMBA != "Timed out!"]

MITISHAMBA_3 <- bind_rows(MITISHAMBA_2)

MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub("\r\n", " ", x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub("\r", " ", x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub("\n", " ", x)
  })
MITISHAMBA_3[] <-
  lapply(MITISHAMBA_3, function(x) {
    gsub("\t", " ", x)
  })

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, MITISHAMBA_3)
