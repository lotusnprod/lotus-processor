#' Title
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
myDirtyF <- function(table) {
  table_tot <- table %>%
    group_by(referenceType) %>%
    count(name = "tot")

  table_y <- table %>%
    filter(validated == "Y") %>%
    group_by(referenceType) %>%
    count(name = "y")

  table_n <- table %>%
    filter(validated == "N") %>%
    group_by(referenceType) %>%
    count(name = "n")

  table_mix <- table %>%
    filter(validated == "Y/N") %>%
    group_by(referenceType) %>%
    count(name = "mix")

  table_full <- left_join(table_tot, table_y)
  table_full <- left_join(table_full, table_n)
  table_full <- left_join(table_full, table_mix) %>%
    replace(is.na(.), 0) %>%
    mutate(ratio = y / tot)

  return(table_full)
}

#' Title
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
myDirtyC <- function(table) {
  table_tot <- table %>%
    count(name = "tot")

  table_y <- table %>%
    filter(validated == "Y") %>%
    count(name = "y")

  table_n <- table %>%
    filter(validated == "N") %>%
    count(name = "n")

  table_mix <- table %>%
    filter(validated == "Y/N") %>%
    count(name = "mix")

  table_full <-
    bind_cols(table_tot, table_y, table_n, table_mix) %>%
    replace(is.na(.), 0) %>%
    mutate(ratio = y / tot)

  return(table_full)
}

#' Title
#'
#' @param table
#' @param title
#' @param yaxismax
#'
#' @return
#' @export
#'
#' @examples
myDirtyP <- function(table, title, yaxismax) {
  fig <-
    plot_ly(
      data = table,
      x = ~referenceType,
      y = ~y,
      type = "bar",
      name = "correct",
      color = I("green")
    ) %>%
    add_trace(
      y = ~n,
      name = "uncorrect",
      color = I("red"),
      text = ~ round(x = ratio, digits = 2),
      textposition = "outside",
      textfont = list(color = I("black"), size = 20)
    ) %>%
    layout(
      title = title,
      yaxis = list(
        title = "Count",
        range = c(0, yaxismax)
      ),
      barmode = "stack"
    )
  return(fig)
}

#' Title
#'
#' @param table
#' @param title
#' @param yaxismax
#'
#' @return
#' @export
#'
#' @examples
myDirtyQ <- function(table, title, yaxismax) {
  fig <-
    plot_ly(
      data = table,
      y = ~y,
      type = "bar",
      name = "correct",
      color = I("green")
    ) %>%
    add_trace(
      y = ~n,
      name = "uncorrect",
      color = I("red"),
      text = ~ round(x = ratio, digits = 2),
      textposition = "outside",
      textfont = list(color = I("black"), size = 20)
    ) %>%
    layout(
      title = title,
      yaxis = list(
        title = "Count",
        range = c(0, yaxismax)
      ),
      barmode = "stack"
    )
  return(fig)
}
