# title: "TipDB cleaneR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/parallel.R")
source("functions/standardizing.R")

library(dplyr)
library(pbmcapply)
library(readr)
library(splitstackshape)
library(stringr)
library(tidyr)
library(jsonlite)

# get paths
database <- databases$get("tipdb")

## files
ids <- read_delim(
  file = database$sourceFiles$mol,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_names = "id"
) %>%
  mutate_all(as.character) %>%
  filter(grepl(pattern = "*.mol", x = id))

X <- ids$id

getTipInChI <- function(i) {
  path <-
    paste(
      "'/Users/rutza/GitLab/opennaturalproductsdb",
      gsub(
        pattern = "..",
        replacement = "",
        x = X[i],
        fixed = TRUE
      ),
      "'",
      sep = ""
    )
  
  path
  
  file <- str_extract(string = path, pattern = "TIP[0-9]{6}.mol")
  
  command <-
    paste("/usr/local/bin/obabel -imol ", path, " -as   -oinchi")
  
  test <- data.frame(file,
            system(command = command,
                   intern = TRUE))
  return(test)
}

i <- 1:length(X)

tipInchiList <- invisible(
  pbmclapply(
    FUN = getTipInChI,
    X = i,
    mc.silent = FALSE,
    mc.cores = numCores,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE,
    ignore.interactive = TRUE
  )
)

tipInchiDf <- bind_rows(tipInchiList)

colnames(tipInchiDf)[2] <- "inchi"
