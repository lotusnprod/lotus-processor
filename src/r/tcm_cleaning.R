#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tcm_cleaning <- function(x) {
  data <- x

  data$newbiologicalsource <-
    gsub(" bulbus", "", data$latin, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" caulis", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" caulis et folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" corolla", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" cortex", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" exocarpium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" exocarpium rubrum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" flos", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium et cacumen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" folium et caulis", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus germinatus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus immaturus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus retinervus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" fructus rotundus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" herba", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" lignum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" medulla", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" pericarpum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" petiolus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" pollen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radicis cortex", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix et rhizoma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" radix preparata", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramulus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramulus cum uncus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" ramus et folium", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma alba", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" rhizoma et radix", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" semen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" semen germinatum", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" spica", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" stamen", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" stigma", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" thallus", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub(" et", "", data$newbiologicalsource, fixed = TRUE)

  # not tuber
  data_final <- data %>%
    mutate(latin = newbiologicalsource) %>%
    select(latin, common, biologicalsource)

  return(data_final)
}
