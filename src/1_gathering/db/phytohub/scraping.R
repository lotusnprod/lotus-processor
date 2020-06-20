#title: "PHYTOHUB scrapeR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/parallel.R")

library(dplyr)
library(pbmcapply)
library(parallel)
library(data.table)
library(splitstackshape) # provides cSplit
library(stringr) # provides str_pad
library(rvest)  # provides read_html

# get paths
database <- databases$get("phytohub")

url <- 'http://phytohub.eu/entries/PHUB'

X <- (1:1975)

getphytohub <- function(X)
{
  tryCatch({
    cd_id <- str_pad (X, 6, pad = "0")
    url_id <- paste(url, cd_id)
    url_id <- gsub("\\s", "", url_id)
    sample <- read_html(url_id)
    scrape1 <-
      html_nodes(sample, xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[2]") %>%
      html_text()
    scrape2 <-
      html_nodes(sample, xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[10]/pre/small") %>%
      html_text()
    scrape3 <-
      html_nodes(sample, xpath = "//*[@id=\"fs\"]/div/div/div/table") %>%
      html_table()
    scrape4 <- scrape3[[1]]
    scrape5 <-
      html_nodes(sample, xpath = "/html/body/main/div[2]/section[1]/div/div[2]/dl/dd[11]/pre") %>%
      html_text()
    
    df <- cbind(scrape1, scrape2, scrape5, scrape4)
    final_df <- data.frame(df)
    return(final_df)
  },
  error = function(e) {
    "Timed out!"
  })
}

PHYTOHUB <- invisible(
  pbmclapply(
    FUN = getphytohub,
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

PHYTOHUB_2 <- PHYTOHUB[PHYTOHUB != "Timed out!"]

PHYTOHUB_3 <- bind_rows(PHYTOHUB_2)

PHYTOHUB_4 <- PHYTOHUB_3 %>%
  select(
    name = scrape1,
    inchi = scrape2,
    smiles = scrape5,
    biologicalsource = Name,
    name_precursor = Precursor,
    biologicalsource_precursor = Food.Source
  ) %>%
  filter(is.na(name_precursor)) %>%
  select(name,
         inchi,
         smiles,
         biologicalsource)

url <- 'http://phytohub.eu/entry_food_sources/'

X <- (1:2975)

getphytohubref <- function(X)
{
  tryCatch({
    cd_id <- X
    url_id <- paste(url, cd_id)
    url_id <- gsub("\\s", "", url_id)
    sample <- read_html(url_id)
    scrape1 <-
      html_nodes(sample, xpath = "/html/body/main/div/h1") %>%
      html_text()
    scrape2 <-
      html_nodes(sample, xpath = "/html/body/main/blockquote") %>%
      html_text()
    
    df <- cbind(scrape1, scrape2)
    
    final_df <- as.data.frame(df)
    
    return(final_df)
  },
  error = function(e) {
    "Timed out!"
  })
}

PHYTOHUB_REF <- invisible(
  pbmclapply(
    FUN = getphytohubref,
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

PHYTOHUB_5 <- PHYTOHUB_REF[PHYTOHUB_REF != "Timed out!"]

PHYTOHUB_6 <- bind_rows(PHYTOHUB_5)

colnames(PHYTOHUB_6) <- c("pair", "reference")

PHYTOHUB_7 <- PHYTOHUB_6 %>%
  mutate(joining_col = gsub("Publications for ", "", pair)) %>%
  mutate(joining_col = gsub("being present in ", "", joining_col))

PHYTOHUB_8 <- PHYTOHUB_4 %>%
  mutate(joining_col = paste(name, biologicalsource, sep = " "))

PHYTOHUB_9 <- full_join(PHYTOHUB_7, PHYTOHUB_8) %>%
  select(name, inchi, smiles, biologicalsource, reference) %>%
  filter(!is.na(name) | !is.na(inchi) | !is.na(smiles))

PHYTOHUB_9$reference <- y_as_na(PHYTOHUB_9$reference, "")

# exporting
database$writeFile(database$sourceFiles$tsv, PHYTOHUB_9)
