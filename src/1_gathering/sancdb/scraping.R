#title: "SANCDB scrapeR"

#loading
##functions
source("../../functions.R")

outpath <- "0_initial_files/SANCDB_scraped.tsv.zip"

url <- 'https://sancdb.rubi.ru.ac.za/compounds/'

X <- (1:1000)

getsanc <- function(X)
{
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
    mc.allow.recursive = TRUE
  )
) %>%
  as.data.table() %>%
  t() %>%
  cSplit("V1", "\n")

#exporting
write.table(
  x = SANCDB,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
