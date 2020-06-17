library(dplyr)

standardizing_original <- function(data_selected,
                                   db,
                                   structure_field)
{
  data_selected[setdiff(c("name",
                          "biologicalsource",
                          "reference"),
                        names(data_selected))] <- NA
  
  data_standard <- data.frame(data_selected) %>%
    mutate(database = db) %>%
    select(database,
           name,
           all_of(structure_field),
           biologicalsource,
           reference) %>%
    distinct_at(vars(all_of(structure_field),
                     biologicalsource),
                .keep_all = TRUE)
  
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\r\n", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\r", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\n", " ", x))
  data_standard[] <-
    lapply(data_standard, function(x)
      gsub("\t", " ", x))
  
  return(data_standard)
}
