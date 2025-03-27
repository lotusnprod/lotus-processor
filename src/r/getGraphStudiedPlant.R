#' Title
#'
#' @param plant
#'
#' @return
#' @export
#'
#' @examples
getGraphStudiedPlant <- function(plant) {
  try({
    mostplant <-
      as.character(inhouseDb_most_plant[
        inhouseDb_most_plant$organism_name == plant,
        1
      ])
    inhouseDb_most_organism_2plot <-
      inhouseDbMeta %>%
      filter(organism_name == mostplant) %>%
      distinct(
        structure_inchikey_2D,
        organism_name,
        database,
        .keep_all = TRUE
      ) %>%
      group_by(organism_name, database) %>%
      count(structure_inchikey_2D) %>%
      ungroup()
    inhouseDb_most_organism_2plot_wide <-
      inhouseDb_most_organism_2plot %>%
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
    inhouseDb_most_organism_2plot_wide <-
      left_join(
        inhouseDb_most_organism_2plot_wide,
        chemo
      ) %>%
      distinct(structure_inchikey_2D, organism_name, .keep_all = TRUE)
    dbnumostorganism <- as.numeric(nrow(
      inhouseDbMeta %>%
        filter(organism_name == mostplant) %>%
        distinct(database)
    ))
    mostsuperclasses <- inhouseDb_most_organism_2plot_wide %>%
      filter(!is.na(structure_taxonomy_npclassifier_01pathway)) %>%
      count(structure_taxonomy_npclassifier_01pathway) %>%
      arrange(desc(n)) %>%
      head(7)
    upset(
      inhouseDb_most_organism_2plot_wide,
      nsets = 10,
      query.legend = "top",
      queries = list(
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1],
              mostsuperclasses[6, 1],
              mostsuperclasses[7, 1]
            )
          ),
          active = TRUE,
          color = "#fdbf6f",
          query.name = mostsuperclasses[7, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1],
              mostsuperclasses[6, 1]
            )
          ),
          active = TRUE,
          color = "#e31a1c",
          query.name = mostsuperclasses[6, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1]
            )
          ),
          active = TRUE,
          color = "#fb9a99",
          query.name = mostsuperclasses[5, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1]
            )
          ),
          active = TRUE,
          color = "#33a02c",
          query.name = mostsuperclasses[4, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1]
            )
          ),
          active = TRUE,
          color = "#b2df8a",
          query.name = mostsuperclasses[3, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1]
            )
          ),
          active = TRUE,
          color = "#1f78b4",
          query.name = mostsuperclasses[2, 1]
        ),
        list(
          query = elements,
          params = list(
            "structure_taxonomy_npclassifier_01pathway",
            mostsuperclasses[1, 1]
          ),
          active = TRUE,
          color = "#a6cee3",
          query.name = mostsuperclasses[1, 1]
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
      set_size.show = TRUE,
      set_size.scale_max = 400
    )
  })
}
