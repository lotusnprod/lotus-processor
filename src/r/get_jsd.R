library(dplyr)
# library(future)
# library(future.apply)
library(philentropy)
# library(progressr)

# source(file = "r/progressr.R")

#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
get_jsd <- function(X) {
  # p <- progressr::progressor(along = xs)
  # future.apply::future_lapply(
  #   future.seed = TRUE,
  #   X = xs,
  #   FUN = function(x) {
  #     p(sprintf("x=%g", as.numeric(x))) ## little hack
  table_real <- table_counted |>
    dplyr::filter(!is.na(!!as.name(chem_level))) %>%
    dplyr::distinct(!!as.name(bio_level), !!as.name(chem_level), !!as.name(paste(chem_level, bio_level,
      sep =
        "_"
    ))) %>%
    dplyr::filter(!!as.name(chem_level) == Y[X]) %>%
    dplyr::right_join(table %>% dplyr::distinct(!!as.name(bio_level))) %>%
    dplyr::distinct(!!as.name(bio_level), !!as.name(paste(chem_level, bio_level,
      sep =
        "_"
    ))) %>%
    dplyr::filter(!is.na(!!as.name(bio_level))) %>%
    dplyr::mutate_all(~ replace(., is.na(.), 0))

  table_random <- table |>
    dplyr::filter(!is.na(!!as.name(bio_level))) |>
    dplyr::distinct(!!as.name(bio_level)) |>
    cbind(1) ## whatever numbre since distribution

  P <- table_real[[2]]
  Q <- table_random[[2]]

  x <- rbind(
    P,
    Q
  )

  jsd <- philentropy::JSD(x = x, est.prob = "empirical")

  result <- data.frame(
    Y[X],
    jsd
  ) |>
    rename(!!as.name(chem_level) := Y.X., !!as.name(bio_level) := jsd)

  return(result)
}
#   )
# }
