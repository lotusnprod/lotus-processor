#' Title
#'
#' @param table
#' @param group
#' @param variable
#'
#' @return
#' @export
#'
#' @examples
prepare_upset <- function(table, group, variable) {
  table_2plot <- table %>%
    distinct(!!!syms(variable), !!!syms(group)) %>%
    group_by(!!!syms(group)) %>%
    count(!!!syms(variable)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = !!as.name(group),
      values_from = n
    ) %>%
    mutate_at(
      .vars = (length(group) + length(variable)):ncol(.),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = (length(group) + length(variable)):ncol(.),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    distinct(!!!syms(variable), .keep_all = TRUE) %>%
    data.frame()

  return(table_2plot)
}

#' Title
#'
#' @param table
#' @param group
#' @param variable
#'
#' @return
#' @export
#'
#' @examples
prepare_upset_complex <- function(table, group, variable) {
  table_2plot <- table %>%
    distinct(!!!syms(variable), !!!syms(group)) %>%
    group_by(!!!syms(group)) %>%
    count(!!!syms(variable)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = group[1],
      values_from = n
    ) %>%
    mutate_at(
      .vars = (length(group) + length(variable)):ncol(.),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = (length(group) + length(variable)):ncol(.),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    distinct(!!as.name(group[2]), .keep_all = TRUE) %>%
    data.frame()

  return(table_2plot)
}
