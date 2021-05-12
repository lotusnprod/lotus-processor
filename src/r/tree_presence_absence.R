#' Title
#'
#' @param table
#' @param level
#'
#' @return
#' @export
#'
#' @examples
tree_presence_absence <- function(table, level) {
  q <-
    ggtree(
      tr = tr_temp,
      layout = "circular",
      size = 5,
      aes(color = switch(level,
        "structure_inchikey" = structure_inchikey,
        "structure_taxonomy_npclassifier_03class" = structure_taxonomy_npclassifier_03class,
        "structure_taxonomy_npclassifier_02superclass" = structure_taxonomy_npclassifier_02superclass,
        "structure_taxonomy_npclassifier_01pathway" = structure_taxonomy_npclassifier_01pathway
      ))
    ) %<+%
    table +
    scale_color_manual(
      values = c("#336A2D", "#861450"),
      na.value = "grey"
    )
  q <- q %<+%
    info_temp +
    new_scale_color() +
    geom_tiplab(
      aes(color = Kingdom),
      align = TRUE,
      size = rel(5),
      offset = rel(1)
    ) +
    scale_color_manual(
      values = strsplit(x = paired, split = " "),
      na.value = "grey"
    ) +
    theme(legend.position = "none")
}
