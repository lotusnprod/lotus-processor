# title: "ALKAMID scrapeR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

url <- 'http://alkamid.ugent.be/molecule.php?ID='

X <- (1:439)

getalkamid <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "&typegroup=genusgroup")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node("body") %>%
      html_node(xpath = "/html/body/div[2]/div/div[2]/div/div[2]/div/div[2]/div/div/div/div/div/div/div") %>%
      html_text()
  },
  error = function(e) {
    NA
  })
}

getalkamid_ref <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, "&typegroup=genusgroup")
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node("body") %>%
      html_node(xpath = "/html/body/div[2]/div/div[2]/div/div[2]/div/div[2]/div/div/div/div/div/div/div/div[8]/div[2]/table") %>%
      html_table() %>%
      mutate_all(as.character)
  },
  error = function(e) {
    NA
  })
}

ALKAMID <- invisible(
  pbmclapply(
    FUN = getalkamid,
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

ALKAMID_REF <- invisible(
  pbmclapply(
    FUN = getalkamid_ref,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

ALKAMID_REF_2 <- ALKAMID_REF[!is.na(ALKAMID_REF)]

ALKAMID_REF_3 <- bind_rows(ALKAMID_REF_2, .id = "entry_id")

# exporting
write.table(
  x = ALKAMID,
  file = gzfile(
    description = pathDataExternalDbSourceAlkamidOriginal,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)

write.table(
  x = ALKAMID_REF_3,
  file = gzfile(
    description  = pathDataExternalDbSourceAlkamidRef,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
