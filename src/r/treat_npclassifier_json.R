#' Title
#'
#' @return
#' @export
#'
#' @examples
treat_npclassifier_json <- function() {
  taxonomy_classes <- taxonomy$Class %>%
    rbind()
  rownames(taxonomy_classes) <- "id_class"
  taxonomy_classes <- taxonomy_classes %>%
    t() %>%
    data.frame() %>%
    mutate(
      class = rownames(.),
      id_class = as.numeric(id_class)
    )

  taxonomy_superclasses <- taxonomy$Superclass %>%
    rbind()
  rownames(taxonomy_superclasses) <- "id_superclass"
  taxonomy_superclasses <- taxonomy_superclasses %>%
    t() %>%
    data.frame() %>%
    mutate(
      superclass = rownames(.),
      id_superclass = as.numeric(id_superclass)
    )

  taxonomy_pathways <- taxonomy$Pathway %>%
    rbind()
  rownames(taxonomy_pathways) <- "id_pathway"
  taxonomy_pathways <- taxonomy_pathways %>%
    t() %>%
    data.frame() %>%
    mutate(
      pathway = rownames(.),
      id_pathway = as.numeric(id_pathway)
    )

  taxonomy_hierarchy_class <- taxonomy$Class_hierarchy

  id_pathway <- list()
  id_superclass <- list()
  id_class <- list()

  for (i in seq_len(length(taxonomy_hierarchy_class))) {
    id_pathway[[i]] <- taxonomy_hierarchy_class[[i]]$Pathway
    id_superclass[[i]] <- taxonomy_hierarchy_class[[i]]$Superclass
    id_class[[i]] <- names(taxonomy_hierarchy_class[i])
  }
  zu <- cbind(id_pathway, id_superclass, id_class) %>%
    data.frame() %>%
    mutate(id_class = as.numeric(id_class)) %>%
    unnest(id_superclass) %>%
    unnest(id_pathway)

  ## No idea why would this be needed... class already has everything?
  id_pathway_2 <- list()
  id_superclass <- list()

  taxonomy_hierarchy_superclass <- taxonomy$Super_hierarchy

  for (i in seq_len(length(taxonomy_hierarchy_superclass))) {
    id_pathway_2[[i]] <- taxonomy_hierarchy_superclass[[i]]$Pathway
    id_superclass[[i]] <- names(taxonomy_hierarchy_superclass[i])
  }
  zu_2 <- cbind(id_pathway_2, id_superclass) %>%
    data.frame() %>%
    mutate(id_superclass = as.numeric(id_superclass)) %>%
    unnest(id_pathway_2)

  taxonomy_semicleaned <- full_join(zu, taxonomy_classes) %>%
    full_join(., taxonomy_superclasses) %>%
    full_join(., taxonomy_pathways) %>%
    distinct(class, superclass, pathway)

  return(taxonomy_semicleaned)
}
