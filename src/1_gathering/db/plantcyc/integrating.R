#title: "PLANTCYC compileR"

# loading paths
source("paths.R")
source("functions/helpers.R")
source("functions/standardizing.R")

library(dplyr)
library(readr)
library(splitstackshape)
library(tidyr)

# get paths
database <- databases$get("plantcyc")

# defining plantcycompiling function
plantcycompiling <- function(x)
{
  path <- paste(x,
                "/",
                file,
                sep = "")
  
  data <- read.delim(file = path)
  
  data$X. <- as.character(data$X.)
  
  data$X. <- sub(" - ", " ", data$X.)
  
  data_2 <- data %>%
    cSplit(splitCols = "X.",
           sep = " ") %>%
    select(1:2) %>%
    tibble()
  
  colnames(data_2)[1] <- "A"
  colnames(data_2)[2] <- "B"
  
  data_2$A <- as.character(data_2$A)
  
  data_3 <- data_2 %>%
    filter(A == "INCHI" |
             A == "UNIQUE-ID")
  
  data_4 <- data_3 %>%
    filter(A == "UNIQUE-ID" &
             lead(A, n = 1) == "INCHI" |
             A == "INCHI" &
             lag(A, n = 1) == "UNIQUE-ID")
  
  data_5 <- data_4 %>%
    pivot_wider(names_from = A,
                values_from = B) %>%
    unnest()
  
  data$X. <- as.character(data$X.)
  
  test <- data %>%
    filter(., grepl(pattern = "# Organism: ",
                    x = X.))
  
  test_2 <-
    data.frame(gsub(pattern = "# Organism: ",
                    replacement = "",
                    test)) %>% cSplit(splitCols = 1,
                                      sep = " ")
  
  colnames(test_2)[1] <- "A"
  
  test_2$A <- as.character(test_2$A)
  
  colnames(test_2)[2] <- "B"
  
  test_2$B <- as.character(test_2$B)
  
  biologicalsource <- paste(test_2$A,
                            test_2$B)
  
  colnames(data_5)[1] <- "uniqueid"
  
  colnames(data_5)[2] <- "inchi"
  
  data_5$biologicalsource <- biologicalsource
  
  data_standard <- data_5
  
  biologicalsource_2 <- gsub(pattern = " ",
                             replacement = "_",
                             x = biologicalsource)
  
  outpath <- file.path(pathDataExternalDbSourcePlantcyc,
                       paste(biologicalsource_2,
                             ".tsv.zip",
                             sep = ""))
  
  write.table(
    x = data_standard,
    file = gzfile(
      description = outpath,
      compression = 9,
      encoding = "UTF-8"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    fileEncoding = "UTF-8"
  )
}


## paths
dirnames <-
  data.frame(list.dirs(pathDataExternalDbSourcePlantcycDir))

dirnames[, 1] <- as.character(dirnames[, 1])

dirnames <- dirnames[seq(4, nrow(dirnames), 6), ]

file <- "compounds.dat"

for (i in dirnames) {
  plantcycompiling(x = i)
}
