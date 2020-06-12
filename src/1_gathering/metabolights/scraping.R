#title: "Metabolights studies scrapeR"

# setting working directory
setwd("~/GitLab/opennaturalproductsdb/src/")

# loading paths
source("paths.R")

# loading functions
source("functions.R")

data <- xmlParse(pathDataExternalDbSourceMetabolightsStudies)

xml_data <- xmlToList(data)

n <- length(xml_data[["entries"]][[1]])
df <-
  structure(xml_data[["entries"]], row.names = c(NA, -n), class = "data.frame")

entries <- as.data.frame(unlist(df[3,]))

entries_filtered <- entries %>%
  filter(grepl("MTBLC", `unlist(df[3, ])`)) %>%
  distinct(.)

url <-
  'https://www.ebi.ac.uk/metabolights/webservice/beta/compound/'

X <- entries_filtered$`unlist(df[3, ])`

getmetabolights <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id, sep = "")
    destfile = paste(pathDataExternalDbSourceMetabolightsStudiesScrapedDir,
                     cd_id,
                     ".json",
                     sep = "")
    text <- read_html(url_id) %>%
      html_text()
  },
  error = function(e) {
    
  })
  write(text,
        file = destfile)
}

pbmclapply(
  FUN = getmetabolights,
  X = X,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)
