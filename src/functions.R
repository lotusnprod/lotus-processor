#######################################################
####################   Functions   ####################
#######################################################

library(data.table)
library(dplyr)
library(parallel)
library(pbmcapply)
library(readr)
library(rvest)
library(stringr)
library(tidyr)

#######################################################
#######################################################

db_loader <- function(path_to_db) {
  db <- read_delim(
    file = gzfile(path_to_db),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  return(db)
}

#######################################################
#######################################################

shift <- function(x, n) {
  c(x[-(seq(n))], rep(NA, n))
}

#######################################################
#######################################################

split_data_table <-
  function(x, no_rows_per_frame, text, path_to_store) {
    split_vec <- seq(1, nrow(x), no_rows_per_frame)
    
    for (split_cut in split_vec) {
      sample <- x[split_cut:(split_cut + (no_rows_per_frame - 1))]
      write.table(
        sample,
        paste(path_to_store,
              text,
              as.integer(split_cut + (
                no_rows_per_frame - 1
              )),
              ".tsv",
              sep = ""),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        fileEncoding = "UTF-8"
      )
    }
  }

#######################################################
#######################################################

y_as_na <- function(x, y)
{
  if ("factor" %in% class(x))
    x <- as.character(x) ## since ifelse wont work with factors
  ifelse(test = as.character(x) != y,
         yes = x,
         no = NA)
}
