# title: "NANPDB scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

groundhog.library(dplyr, date = groundhog.day)
groundhog.library(pbmcapply, date = groundhog.day)
library(parallel)
groundhog.library(data.table, date = groundhog.day)
groundhog.library(splitstackshape, date = groundhog.day) # provides cSplit
groundhog.library(rvest, date = groundhog.day) # provides read_html

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
        html_node("body") %>%
        html_node(xpath = "/html/body/div[3]/div/div[2]") %>%
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
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
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
database$writeFile(database$sourceFiles$tsv, NANPDB_2)
