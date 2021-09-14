# title: "ALKAMID scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

library(data.table)
library(dplyr)
library(parallel)
library(pbmcapply)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

# get paths
database <- databases$get("alkamid")

url <- "http://alkamid.ugent.be/molecule.php?ID="

X <- (1:439)

getalkamid <- function(X) {
  cd_id <- X
  url_id <- paste(url, cd_id, "&typegroup=genusgroup")
  url_id <- gsub("\\s", "", url_id)
  html <- read_html(url_id) %>% html_element("body")
  data <-
    html %>%
    html_element(xpath = "/html/body/div[2]/div/div[2]/div/div[2]/div/div[2]/div/div/div/div/div/div/div") %>%
    html_text()
  ref <- html %>%
    html_element(xpath = "/html/body/div[2]/div/div[2]/div/div[2]/div/div[2]/div/div/div/div/div/div/div/div[8]/div[2]/table") %>%
    html_table() %>%
    mutate_all(as.character)
  list(data = data, ref = ref)
}

extracted_elements <- invisible(
  pbmclapply(
    FUN = getalkamid,
    X = X,
    mc.silent = FALSE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

ALKAMID <- lapply(extracted_elements, function(x) {
  x$data
}) %>%
  as.data.table() %>%
  t() %>%
  cSplit("V1", "\n")

ALKAMID_REF <- lapply(extracted_elements, function(x) {
  x$ref
})
ALKAMID_REF_2 <- ALKAMID_REF[!is.na(ALKAMID_REF)]
ALKAMID_REF_3 <- bind_rows(ALKAMID_REF_2, .id = "entry_id")

# exporting
ifelse(
  test = !dir.exists(dirname(database$sourceFiles$tsv)),
  yes = dir.create(dirname(database$sourceFiles$tsv)),
  no = paste(dirname(database$sourceFiles$tsv), "exists")
)

database$writeFile(database$sourceFiles$tsv, ALKAMID)
database$writeFile(database$sourceFiles$tsvRef, ALKAMID_REF_3)
