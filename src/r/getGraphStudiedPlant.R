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
      as.character(inhouseDb_most_plant[inhouseDb_most_plant$organismCleaned == plant, 1])
    inhouseDb_most_organism_2plot <-
      inhouseDbMeta %>%
      filter(organismCleaned == mostplant) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        database,
        .keep_all = TRUE
      ) %>%
      group_by(organismCleaned, database) %>%
      count(structureCleanedInchikey2D) %>%
      ungroup()
    inhouseDb_most_organism_2plot_wide <-
      inhouseDb_most_organism_2plot %>%
      pivot_wider(
        names_from = database,
        values_from = n
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = is.na(.),
          values = 0
        )
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
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
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        .keep_all = TRUE
      )
    dbnumostorganism <- as.numeric(nrow(
      inhouseDbMeta %>%
        filter(organismCleaned == mostplant) %>%
        distinct(database)
    ))
    mostsuperclasses <- inhouseDb_most_organism_2plot_wide %>%
      filter(!is.na(structureCleaned_classyfire_2superclass)) %>%
      count(structureCleaned_classyfire_2superclass) %>%
      arrange(desc(n)) %>%
      head(10)
    upset(
      inhouseDb_most_organism_2plot_wide,
      nsets = 10,
      query.legend = "top",
      queries = list(
        list(
          query = elements,
          params = list(
            "structureCleaned_classyfire_2superclass",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1],
              mostsuperclasses[6, 1],
              mostsuperclasses[7, 1],
              mostsuperclasses[8, 1],
              mostsuperclasses[9, 1],
              mostsuperclasses[10, 1]
            )
          ),
          active = TRUE,
          color = "#6a3d9a",
          query.name = mostsuperclasses[10, 1]
        ),
        list(
          query = elements,
          params = list(
            "structureCleaned_classyfire_2superclass",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1],
              mostsuperclasses[6, 1],
              mostsuperclasses[7, 1],
              mostsuperclasses[8, 1],
              mostsuperclasses[9, 1]
            )
          ),
          active = TRUE,
          color = "#cab2d6",
          query.name = mostsuperclasses[9, 1]
        ),
        list(
          query = elements,
          params = list(
            "structureCleaned_classyfire_2superclass",
            c(
              mostsuperclasses[1, 1],
              mostsuperclasses[2, 1],
              mostsuperclasses[3, 1],
              mostsuperclasses[4, 1],
              mostsuperclasses[5, 1],
              mostsuperclasses[6, 1],
              mostsuperclasses[7, 1],
              mostsuperclasses[8, 1]
            )
          ),
          active = TRUE,
          color = "#ff7f00",
          query.name = mostsuperclasses[8, 1]
        ),
        list(
          query = elements,
          params = list(
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
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
            "structureCleaned_classyfire_2superclass",
            c(mostsuperclasses[1, 1])
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
      set_size.show = TRUE
    )
  })
}
