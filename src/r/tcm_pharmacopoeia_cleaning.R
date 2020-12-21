" Title
#"
#' @param x
#'
#' @return
#' @export
#'
#' @examples
tcm_pharmacopoeia_cleaning <- function(x) {
  data <- x

  data$newbiologicalsource <-
    gsub("Bulbus ", "", data$organismTranslated, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Cacumen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Caulis ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Corolla ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Cortex ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Exocarpium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Exocarpium rubrum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Flos ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Folium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Folium ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus germinatus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus immaturus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus retinervus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Fructus rotundus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Herba ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Lignum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Medulla ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Pericarpum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Petiolus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Pollen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radicis cortex ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radix ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Radix et rhizoma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Ramulus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Ramus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Rhizoma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Semen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Semen germinatum ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Spica ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Stamen ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Stigma ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Thallus ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("Uncis ", "", data$newbiologicalsource, fixed = TRUE)

  data$newbiologicalsource <-
    gsub("et ", "", data$newbiologicalsource, fixed = TRUE)

  # not tuber
  data_final <- data %>%
    select(-organismTranslated) %>%
    mutate(organismTranslated = newbiologicalsource) %>%
    select(-newbiologicalsource)

  return(data_final)
}
