#' Title
#'
#' @param subclass
#'
#' @return
#' @export
#'
#' @examples
getGraphChemicalClass <- function(subclass) {
  inhouseDb_most_chemical_class_2plot <-
    inhouseDbMeta %>%
    filter(structure_taxonomy_npclassifier_03class == subclass) %>%
    distinct(
      structure_inchikey_2D,
      organism_name,
      database,
      .keep_all = TRUE
    ) %>%
    group_by(organism_taxonomy_06family, database) %>%
    count(structure_inchikey_2D) %>%
    ungroup()

  inhouseDb_most_chemical_class_2plot_wide <-
    inhouseDb_most_chemical_class_2plot %>%
    pivot_wider(
      names_from = database,
      values_from = n
    ) %>%
    mutate_at(
      .vars = 3:ncol(.),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = 3:ncol(.),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    data.frame()

  dbnumostchemicalclass <- as.numeric(nrow(
    inhouseDbMeta %>%
      filter(structure_taxonomy_npclassifier_03class == subclass) %>%
      distinct(database)
  ))

  mostfamilies <-
    inhouseDb_most_chemical_class_2plot_wide %>%
    filter(!is.na(organism_taxonomy_06family)) %>%
    count(organism_taxonomy_06family) %>%
    arrange(desc(n)) %>%
    head(10)

  upset(
    inhouseDb_most_chemical_class_2plot_wide,
    nsets = 10,
    query.legend = "top",
    queries = list(
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1],
            mostfamilies[6, 1],
            mostfamilies[7, 1],
            mostfamilies[8, 1],
            mostfamilies[9, 1],
            mostfamilies[10, 1]
          )
        ),
        active = TRUE,
        color = "#6a3d9a",
        query.name = mostfamilies[10, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1],
            mostfamilies[6, 1],
            mostfamilies[7, 1],
            mostfamilies[8, 1],
            mostfamilies[9, 1]
          )
        ),
        active = TRUE,
        color = "#cab2d6",
        query.name = mostfamilies[9, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1],
            mostfamilies[6, 1],
            mostfamilies[7, 1],
            mostfamilies[8, 1]
          )
        ),
        active = TRUE,
        color = "#ff7f00",
        query.name = mostfamilies[8, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1],
            mostfamilies[6, 1],
            mostfamilies[7, 1]
          )
        ),
        active = TRUE,
        color = "#fdbf6f",
        query.name = mostfamilies[7, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1],
            mostfamilies[6, 1]
          )
        ),
        active = TRUE,
        color = "#e31a1c",
        query.name = mostfamilies[6, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1],
            mostfamilies[5, 1]
          )
        ),
        active = TRUE,
        color = "#fb9a99",
        query.name = mostfamilies[5, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1],
            mostfamilies[4, 1]
          )
        ),
        active = TRUE,
        color = "#33a02c",
        query.name = mostfamilies[4, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1],
            mostfamilies[3, 1]
          )
        ),
        active = TRUE,
        color = "#b2df8a",
        query.name = mostfamilies[3, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          c(
            mostfamilies[1, 1],
            mostfamilies[2, 1]
          )
        ),
        active = TRUE,
        color = "#1f78b4",
        query.name = mostfamilies[2, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_06family",
          mostfamilies[1, 1]
        ),
        active = TRUE,
        color = "#a6cee3",
        query.name = mostfamilies[1, 1]
      )
    ),
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 20,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per database",
    set_size.show = TRUE
  )
}
