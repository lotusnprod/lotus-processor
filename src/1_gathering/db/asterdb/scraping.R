# title: "ASTERDB scrapeR"

# loading paths
source("paths.R")
source("r/parallel.R")

library(tidyverse)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest) # provides read_html

# get paths
database <- databases$get("asterdb")

url <- "http://143.107.203.160:8080/search/molecule?id="

X <- (1:20000)

getaster <- function(X) {
  tryCatch(
    {
      cd_id <- X
      url_id <- paste0(url, cd_id)
      html <- read_html(url_id) %>% html_node("body")
      data <- html %>%
        xml_child(2) %>%
        xml_child(3) %>%
        xml_child(1) %>%
        html_table(fill = TRUE) %>%
        pivot_wider(names_from = X1, values_from = X2)
      ref <- html %>%
        xml_child(2) %>%
        xml_child(4) %>%
        xml_child(2) %>%
        html_table(fill = TRUE)
      return(cbind(data, ref))
    },
    error = function(e) {
      NA
    }
  )
}

extracted_elements <- invisible(
  pbmclapply(
    FUN = getaster,
    X = X,
    mc.silent = FALSE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

ASTERDB <- bind_rows(extracted_elements[!is.na(extracted_elements)])

## # exporting
database$writeFile(database$sourceFiles$tsv, ALKAMID)
database$writeFile(database$sourceFiles$tsvRef, ALKAMID_REF_3)
