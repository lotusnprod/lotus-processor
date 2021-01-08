# title: "SANCDB scrapeR"

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
database <- databases$get("sancdb")

url <- "https://sancdb.rubi.ru.ac.za/compounds/"

X <- (1:1000)

getsanc <- function(X) {
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "/")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node("body") %>%
      html_node("div#wrap") %>%
      html_node("div#content.content") %>%
      html_node("div#pt-main.pt-perspective") %>%
      html_text()
  })
}

SANCDB <- invisible(
  pbmclapply(
    FUN = getsanc,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
) %>%
  as.data.table() %>%
  t() %>%
  cSplit("V1", "\n")

# exporting
database$writeFile(database$sourceFiles$tsv, SANCDB)