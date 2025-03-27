#' Title
#'
#' @param data
#' @param biological_level
#' @param chemical_level
#' @param biological_filter_value
#' @param biological_filter_level
#' @param chemical_filter_value
#' @param chemical_filter_level
#' @param palette
#'
#' @return
#' @export
#'
#' @examples
draw_chord <-
  function(
    data,
    biological_level,
    chemical_level,
    biological_filter_value = NULL,
    biological_filter_level = NULL,
    chemical_filter_value = NULL,
    chemical_filter_level = NULL,
    palette = paired_palette_med
  ) {
    try({
      table <- data.frame(data)
      table <- table[!is.na(table[, biological_level]), ]
      table <- table[!is.na(table[, chemical_level]), ]
      if (!is.null(biological_filter_value)) {
        table <-
          table[table[, biological_filter_level] %in% biological_filter_value, ]
      }
      if (!is.null(chemical_filter_value)) {
        table <-
          table[table[, chemical_filter_level] %in% chemical_filter_value, ]
      }
      m1 <-
        as.data.table(table(table[, c(biological_level, chemical_level)]))
      m1 <- m1 %>%
        pivot_wider(names_from = 2, values_from = N) %>%
        unnest() %>%
        column_to_rownames(var = biological_level) %>%
        select(order(colSums(-.)))
      m2 <- t(m1) %>%
        as.data.frame() %>%
        select(order(colSums(-.)))
      m2$name <- colnames(m1)
      # colnames(m1) <- paste("chemo", colnames(m1), sep = "_")
      m1$name <- sort(colnames(m2[, seq_len(ncol(m2)) - 1]))
      # colnames(m2)[1:ncol(m2)-1] <- paste("bio", colnames(m2[1:ncol(m2)-1]), sep = "_")
      test_3 <- full_join(m2, m1)
      test_3[is.na(test_3)] <- 0
      rownames(test_3) <- test_3$name
      test_3 <- test_3 %>% select(-name)
      test_4 <- as.matrix(test_3)
      test_5 <- test_4[colnames(test_4), colnames(test_4)]
      chord <- chorddiag(
        data = test_5,
        groupColors = palette,
        groupnamePadding = 10,
        groupThickness = 0.1,
        chordedgeColor = palette,
        groupnameFontsize = 18,
        margin = 260,
        showTooltips = FALSE,
        ticklabelFontsize = 0,
        showTicks = FALSE,
        showZeroTooltips = FALSE
      )
      return(chord)
    })
  }
