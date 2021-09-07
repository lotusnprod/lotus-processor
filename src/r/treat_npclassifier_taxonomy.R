#' Title
#'
#' @return
#' @export
#'
#' @examples
treat_npclassifier_taxonomy <- function() {
  chemical_taxonomy_1 <- chemical_taxonomy_1 %>%
    filter(!is.na(pathway)) %>%
    filter(
      pathway != "Fatty acids" |
        superclass != "Fatty acyls" |
        class != "Halogenated hydrocarbons" |
        grepl(pattern = "Cl|Br|I|F", x = structure_smiles_2D) ## because actually returned instead of NA/null
    ) %>%
    distinct(
      structure_smiles_2D,
      structure_taxonomy_npclassifier_01pathway = pathway,
      structure_taxonomy_npclassifier_02superclass = superclass,
      structure_taxonomy_npclassifier_03class = class
    ) %>%
    group_by(structure_smiles_2D) %>%
    summarize(across(
      c(
        "structure_taxonomy_npclassifier_01pathway",
        "structure_taxonomy_npclassifier_02superclass",
        "structure_taxonomy_npclassifier_03class"
      ),
      ~ gsub(
        pattern = "\\bNA\\b",
        replacement = "",
        x = paste(unique(.x), collapse = "|")
      )
    )) %>%
    ungroup() %>%
    mutate(across(
      everything(),
      ~ y_as_na(x = .x, y = "")
    )) %>%
    distinct()

  return(chemical_taxonomy_1)
}
