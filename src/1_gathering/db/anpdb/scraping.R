# title: "ANPDB scrapeR"

# loading paths
source("paths.R")

library(dplyr)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html
library(xml2)

# get paths
database <- databases$get("anpdb")

url <- "http://african-compounds.org/anpdb/get_compound_card/"

X <- 1:8410

getnanp <- function(X) {
  tryCatch(
    {
      cd_id <- X
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

NANPDB <- invisible(lapply(
  FUN = getnanp,
  X = X
))

NANPDB_2 <- bind_rows(NANPDB[!is.na(NANPDB)])

NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "\r\n",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "\r",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "\n",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "\t",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
    replacement = " ",
    x = x
  )
})
NANPDB_2[] <- lapply(NANPDB_2, function(x) {
  gsub(
    pattern = "  ",
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

database$writeFile(database$sourceFiles$tsv, NANPDB_2)
