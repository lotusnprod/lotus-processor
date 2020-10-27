#######################################################
####################   Functions   ####################
#######################################################

library(tidyverse)

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
  return(db)()
}

#######################################################
#######################################################

split_data_table <-
  function(x, no_rows_per_frame, text, path_to_store) {
    split_vec <- seq(1, nrow(x), no_rows_per_frame)
    
    for (split_cut in split_vec) {
      sample <- x[split_cut:(split_cut + (no_rows_per_frame - 1))]
      
      sample <- sample %>%
        filter(!is.na(sample[, 1]))
      
      write.table(
        x = sample,
        file = paste(
          path_to_store,
          text,
          str_pad(
            string = as.integer(split_cut + (no_rows_per_frame - 1)),
            width = 6,
            pad = "0"
          ),
          ".tsv",
          sep = ""
        ),
        row.names = FALSE,
        quote = TRUE,
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
