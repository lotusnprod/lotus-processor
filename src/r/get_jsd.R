library(dplyr)
library(future)
library(future.apply)
library(philentropy)
library(progressr)

source(file = "r/progressr.R")

#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
get_jsd <- function(x) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    future.seed = TRUE,
    X = xs,
    FUN = function(x) {
      p(sprintf("x=%g", as.numeric(x))) ## little hack
      table_real <- table_counted |>
        filter(!is.na(!!as.name(chem_level))) |>
        distinct(!!as.name(bio_level), !!as.name(chem_level), !!as.name(paste(chem_level, bio_level,
          sep =
            "_"
        ))) |>
        filter(!!as.name(chem_level) == Y[x]) |>
        right_join(table |> distinct(!!as.name(bio_level))) |>
        distinct(!!as.name(bio_level), !!as.name(paste(chem_level, bio_level,
          sep =
            "_"
        ))) |>
        filter(!is.na(!!as.name(bio_level))) |>
        mutate_all(~ replace(., is.na(.), 0))

      table_random <- table |>
        filter(!is.na(!!as.name(bio_level))) |>
        distinct(!!as.name(bio_level)) |>
        cbind(1) ## whatever numbre since distribution

      P <- table_real[[2]]
      Q <- table_random[[2]]

      x <- rbind(
        P,
        Q
      )

      jsd <- JSD(x, est.prob = "empirical")

      result <- data.frame(
        Y[x],
        jsd
      ) |>
        rename(!!as.name(chem_level) := Y.x., !!as.name(bio_level) := jsd)

      return(result)
    }
  )
}
