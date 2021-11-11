# title: "NANPDB scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("nanpdb")

url <- "http://african-compounds.org/nanpdb/get_compound_card/"

X <- (1:8410)

getnanp <- function(X) {
  tryCatch(
    {
      cd_id <- X
      url_id <- paste(url, cd_id, "/")
      url_id <- gsub("\\s", "", url_id)
      df1 <- read_html(url_id) %>%
        html_element("body") %>%
        html_element(xpath = "/html/body/div[3]/div/div[2]") %>%
        xml_child(2) %>%
        xml_child(1) %>%
        html_table()

      df2 <- t(df1)

      colnames(df2) <- df2[1, ]

      df3 <- data.frame(df2) %>%
        filter(rownames(.) == "X2")

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

NANPDB <- invisible(
  pbmclapply(
    FUN = getnanp,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE,
    mc.style = "txt",
    mc.substyle = 1
  )
)

NANPDB_2 <- bind_rows(NANPDB[!is.na(NANPDB)])

NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("\r\n", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("\r", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("\n", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("\t", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub("  ", " ", x)
})

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, NANPDB_2)
