library(stringr)

#' Title
#'
#' @param x
#' @param no_rows_per_frame
#' @param text
#' @param path_to_store
#'
#' @return
#' @export
#'
#' @examples
split_data_table <-
  function(x, no_rows_per_frame, text, path_to_store) {
    split_vec <- seq(1, nrow(x), no_rows_per_frame)

    for (split_cut in split_vec) {
      sample <- x[split_cut:(split_cut + (no_rows_per_frame - 1))]

      sample <- sample %>%
        filter(!is.na(!!as.name(colnames(.))))

      write.table(
        x = sample,
        file = paste0(
          path_to_store,
          "/",
          text,
          str_pad(
            string = as.integer(split_cut + (no_rows_per_frame - 1)),
            width = 6,
            pad = "0"
          ),
          ".tsv"
        ),
        na = "",
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        fileEncoding = "UTF-8"
      )
    }
  }
