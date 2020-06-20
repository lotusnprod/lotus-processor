# title: "TMDB scrapeR"

# loading paths
source("paths.R")
source("functions/parallel.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(rvest)  # provides read_html

# get paths
database <- databases$get("tmdb")

url <- 'http://pcsb.ahau.edu.cn:8080/TCDB/f/browseDetail?id='

X <- (1:1473)

gettmdb <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id)
    url_id <- gsub("\\s", "", url_id)
    sample <- read_html(url_id)
    scrape1 <-
      html_nodes(sample, xpath = "/html/body/div[1]/div/table") %>%
      html_table(., fill = TRUE)
    
    scrape2 <- scrape1[[1]]
    return(scrape2)
  },
  error = function(e) {
    "Timed out!"
  })
}

TMDB <- invisible(
  pbmclapply(
    FUN = gettmdb,
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

TMDB_2 <- TMDB[TMDB != "Timed out!"]

TMDB_3 <- bind_rows(TMDB_2)

TMDB_4 <- TMDB_3 %>%
  filter(!is.na(X1))

# exporting
database$writeFile(database$sourceFiles$tsv, TMDB_4)
