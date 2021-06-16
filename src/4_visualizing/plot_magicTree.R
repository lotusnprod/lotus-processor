source("r/log_debug.R")
log_debug("This script plots the magic tree.")

start <- Sys.time()

library(data.table)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
library(rotl)
library(splitstackshape)
library(tidyverse)
source(file = "paths.R")
source(file = "r/colors.R")
source(file = "r/tree_presence_absence.R")

pairs_metadata <- fread(
  input = file.path(
    pathDataProcessed,
    pathLastFrozen
  ),
  na.strings = ""
)

n_min <- 50

families_restricted <- pairs_metadata %>%
  filter(!is.na(organism_taxonomy_06family)) %>%
  group_by(organism_taxonomy_06family) %>%
  add_count() %>%
  ungroup() %>%
  filter(n >= n_min) %>%
  distinct(organism_taxonomy_06family)

families_matched_restricted <-
  tnrs_match_names(
    names = families_restricted$organism_taxonomy_06family,
    do_approximate_matching = FALSE
  )

ott_in_tree <-
  ott_id(families_matched_restricted)[is_in_tree(ott_id(families_matched_restricted))]

families_restricted <- families_restricted %>%
  filter(organism_taxonomy_06family %in% names(ott_in_tree))

families_matched_restricted <-
  tnrs_match_names(
    names = families_restricted$organism_taxonomy_06family,
    do_approximate_matching = FALSE
  )

families_matched_restricted <- families_matched_restricted %>%
  mutate(key = paste(
    gsub(
      x = unique_name,
      pattern = " ",
      replacement = "_",
      fixed = TRUE
    ),
    paste0("ott", ott_id),
    sep = "_"
  ))

tr_restricted <- tol_induced_subtree(ott_ids = ott_in_tree)

specific_classes <- pairs_metadata %>%
  cSplit(
    splitCols = colnames(.)[.[, grepl(
      pattern = "structure_taxonomy_npclassifier_",
      x = colnames(.)
    )]],
    sep = "|",
    direction = "long"
  ) %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_01pathway) &
      !is.na(structure_taxonomy_npclassifier_02superclass) &
      !is.na(structure_taxonomy_npclassifier_03class)
  ) %>%
  mutate_all(as.character) %>%
  filter(organism_taxonomy_06family %in% families_matched_restricted$unique_name) %>%
  filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
  distinct(organism_taxonomy_06family,
    structure_inchikey,
    .keep_all = TRUE
  )

specific_classes_n <- specific_classes %>%
  group_by(organism_taxonomy_06family) %>%
  count(name = "n") %>%
  ungroup()

specific_classes_m <- specific_classes %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  count(name = "m") %>%
  ungroup()

specific_classes_o <- specific_classes %>%
  group_by(organism_taxonomy_06family) %>%
  distinct(organism_taxonomy_08genus, .keep_all = TRUE) %>%
  count(name = "o") %>%
  ungroup()

specific_classes_p <- specific_classes %>%
  group_by(
    organism_taxonomy_06family,
    structure_taxonomy_npclassifier_03class
  ) %>%
  count(name = "p") %>%
  ungroup()

specific_classes <- specific_classes %>%
  left_join(., specific_classes_n) %>%
  left_join(., specific_classes_m) %>%
  left_join(., specific_classes_o) %>%
  left_join(., specific_classes_p) %>%
  mutate(q = p^2 / (m * n)) %>%
  distinct(
    organism_taxonomy_06family,
    structure_taxonomy_npclassifier_03class,
    q,
    o,
    .keep_all = TRUE
  ) %>%
  group_by(organism_taxonomy_06family) %>%
  filter(q == max(q)) %>%
  ungroup() %>%
  distinct(organism_taxonomy_06family, q, .keep_all = TRUE)

tr_restricted$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_restricted$tip.label
  )

taxonomy <-
  left_join(
    families_restricted,
    pairs_metadata
  ) %>%
  distinct(
    Domain = organism_taxonomy_01domain,
    Kingdom = organism_taxonomy_02kingdom,
    Phylum = organism_taxonomy_03phylum,
    Class = organism_taxonomy_04class,
    Order = organism_taxonomy_05order,
    Family = organism_taxonomy_06family
  )

info <- taxonomy %>%
  select(
    id = Family,
    everything()
  ) %>%
  mutate(Kingdom = fct_reorder(Kingdom, !is.na(Domain))) %>%
  mutate(Phylum = fct_reorder(Phylum, !is.na(Kingdom)))

bar_data <- specific_classes %>%
  select(
    id = organism_taxonomy_06family,
    everything(),
    -organism_taxonomy_01domain,
    -organism_taxonomy_02kingdom,
    -organism_taxonomy_03phylum,
    -organism_taxonomy_05order
  ) %>%
  # mutate(`02superclass` = fct_reorder(`02superclass`, `01kingdom`)) %>%
  # mutate(`03class` = fct_reorder(`03class`, !is.na(`02superclass`)))
  mutate(
    superclass = fct_reorder(
      structure_taxonomy_npclassifier_02superclass,
      structure_taxonomy_npclassifier_01pathway
    )
  ) %>%
  mutate(class = fct_reorder(
    structure_taxonomy_npclassifier_03class,
    !is.na(structure_taxonomy_npclassifier_02superclass)
  ))

bar_data_temp <- bar_data

info_temp <- info %>%
  filter(id %in% bar_data_temp$id)

ott_temp <-
  ott_in_tree[names(ott_in_tree) %in% bar_data_temp$id]

ott_temp <- ott_temp[!duplicated(ott_temp)]

tr_temp <- tol_induced_subtree(ott_ids = ott_temp)

tr_temp$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_temp$tip.label
  )

sitosterol_3D <- pairs_metadata %>%
  mutate(structure_inchikey = ifelse(
    test = grepl(
      pattern = "KZJWDPNRJALLNS-VJSFXXLFSA-N",
      x = structure_inchikey,
      fixed = TRUE
    ),
    yes = "KZJWDPNRJALLNS-VJSFXXLFSA-N",
    no = NA
  )) %>%
  group_by(organism_taxonomy_06family) %>%
  fill(structure_inchikey, .direction = "downup") %>%
  distinct(organism_taxonomy_06family, structure_inchikey) %>%
  left_join(families_matched_restricted,
    .,
    by = c("unique_name" = "organism_taxonomy_06family")
  ) %>%
  mutate(structure_inchikey = ifelse(
    test = !is.na(structure_inchikey),
    yes = structure_inchikey,
    no = "noSitosterol_3D"
  )) %>%
  mutate(key = gsub(
    pattern = "_.*",
    replacement = "",
    x = key
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

sitosterol_2D <- pairs_metadata %>%
  mutate(structure_inchikey = ifelse(
    test = grepl(
      pattern = "KZJWDPNRJALLNS",
      x = structure_inchikey,
      fixed = TRUE
    ),
    yes = "KZJWDPNRJALLNS",
    no = NA
  )) %>%
  group_by(organism_taxonomy_06family) %>%
  fill(structure_inchikey, .direction = "downup") %>%
  distinct(organism_taxonomy_06family, structure_inchikey) %>%
  left_join(families_matched_restricted,
    .,
    by = c("unique_name" = "organism_taxonomy_06family")
  ) %>%
  mutate(structure_inchikey = ifelse(
    test = !is.na(structure_inchikey),
    yes = structure_inchikey,
    no = "noSitosterol_2D"
  )) %>%
  mutate(key = gsub(
    pattern = "_.*",
    replacement = "",
    x = key
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

stigmastanes <- pairs_metadata %>%
  mutate(
    structure_taxonomy_npclassifier_03class = ifelse(
      test = structure_taxonomy_npclassifier_03class == "Stigmastane steroids",
      yes = structure_taxonomy_npclassifier_03class,
      no = NA
    )
  ) %>%
  group_by(organism_taxonomy_06family) %>%
  fill(structure_taxonomy_npclassifier_03class, .direction = "downup") %>%
  distinct(
    organism_taxonomy_06family,
    structure_taxonomy_npclassifier_03class
  ) %>%
  left_join(families_matched_restricted,
    .,
    by = c("unique_name" = "organism_taxonomy_06family")
  ) %>%
  mutate(
    structure_taxonomy_npclassifier_03class = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_03class),
      yes = structure_taxonomy_npclassifier_03class,
      no = "zNo_stigmastane_steroids"
    )
  ) %>%
  mutate(key = gsub(
    pattern = "_.*",
    replacement = "",
    x = key
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

steroids <- pairs_metadata %>%
  mutate(
    structure_taxonomy_npclassifier_02superclass = ifelse(
      test = structure_taxonomy_npclassifier_02superclass == "Steroids",
      yes = structure_taxonomy_npclassifier_02superclass,
      no = NA
    )
  ) %>%
  group_by(organism_taxonomy_06family) %>%
  fill(structure_taxonomy_npclassifier_02superclass, .direction = "downup") %>%
  distinct(
    organism_taxonomy_06family,
    structure_taxonomy_npclassifier_02superclass
  ) %>%
  left_join(families_matched_restricted,
    .,
    by = c("unique_name" = "organism_taxonomy_06family")
  ) %>%
  mutate(
    structure_taxonomy_npclassifier_02superclass = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_02superclass),
      yes = structure_taxonomy_npclassifier_02superclass,
      no = "zNo_steroids"
    )
  ) %>%
  mutate(key = gsub(
    pattern = "_.*",
    replacement = "",
    x = key
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

terpenoids <- pairs_metadata %>%
  mutate(
    structure_taxonomy_npclassifier_01pathway = ifelse(
      test = structure_taxonomy_npclassifier_01pathway == "Terpenoids",
      yes = structure_taxonomy_npclassifier_01pathway,
      no = NA
    )
  ) %>%
  group_by(organism_taxonomy_06family) %>%
  fill(structure_taxonomy_npclassifier_01pathway, .direction = "downup") %>%
  distinct(
    organism_taxonomy_06family,
    structure_taxonomy_npclassifier_01pathway
  ) %>%
  left_join(families_matched_restricted,
    .,
    by = c("unique_name" = "organism_taxonomy_06family")
  ) %>%
  mutate(
    structure_taxonomy_npclassifier_01pathway = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_01pathway),
      yes = structure_taxonomy_npclassifier_01pathway,
      no = "zNo_terpenoids"
    )
  ) %>%
  mutate(key = gsub(
    pattern = "_.*",
    replacement = "",
    x = key
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

tree <- ggtree(tr = tr_temp, layout = "circular")

p <- tree %<+%
  info_temp +
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
  geom_fruit(
    data = bar_data_temp,
    geom = geom_bar,
    mapping =
      aes(
        y = id,
        x = q,
        fill = structure_taxonomy_npclassifier_01pathway
      ),
    offset = rel(0.2),
    pwidth = rel(1.6),
    orientation = "y",
    stat = "identity",
  ) +
  geom_tiplab(
    mapping = aes(
      # label = paste(`02superclass`, `03class`, sep = " - "),
      label = paste(
        structure_taxonomy_npclassifier_02superclass,
        structure_taxonomy_npclassifier_03class,
        sep = " - "
      ),
      x = 33 + 75 * q ## to adapt with the data
    ),
    linetype = "blank",
    linesize = 0,
    align = TRUE,
    size = rel(5),
    offset = rel(12)
  ) +
  scale_fill_discrete(
    name = "Chemical pathway",
    direction = "vertical",
    guide = guide_legend(order = 3)
  ) +
  scale_fill_manual(
    values = strsplit(
      x = paired[4:12],
      split = " "
    ),
    na.value = "grey"
  )

p <- p %<+%
  bar_data_temp +
  geom_tippoint(mapping = aes(color = Kingdom, size = o)) +
  new_scale_fill() +
  scale_fill_discrete(
    name = "Biological kingdom",
    guide = guide_legend(
      order = 1,
      ncol = 2
    )
  ) +
  scale_size_continuous(
    name = "Genera in biological family",
    guide = guide_legend(
      order = 2,
      direction = "horizontal"
    )
  ) +
  theme(
    legend.position = c(0.5, 0.05),
    legend.background = element_rect(fill = NA),
    legend.title = element_text(size = rel(3)),
    legend.text = element_text(size = rel(2)),
  )

q <-
  tree_presence_absence(table = sitosterol_3D, level = "structure_inchikey")

r <-
  tree_presence_absence(table = sitosterol_2D, level = "structure_inchikey")

s <-
  tree_presence_absence(table = stigmastanes, level = "structure_taxonomy_npclassifier_03class")

t <-
  tree_presence_absence(table = steroids, level = "structure_taxonomy_npclassifier_02superclass")

u <-
  tree_presence_absence(table = terpenoids, level = "structure_taxonomy_npclassifier_01pathway")

if (mode == "full") {
  ggsave(
    filename = file.path("../res", "magicTree.pdf"),
    plot = p,
    width = 100,
    height = 100,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path("../res", "tree_sitosterol_3D.pdf"),
    plot = q,
    width = 25,
    height = 25,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path("../res", "tree_sitosterol_2D.pdf"),
    plot = r,
    width = 25,
    height = 25,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path("../res", "tree_stigmastanes.pdf"),
    plot = s,
    width = 25,
    height = 25,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path("../res", "tree_steroids.pdf"),
    plot = t,
    width = 25,
    height = 25,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path("../res", "tree_terpenoids.pdf"),
    plot = u,
    width = 25,
    height = 25,
    units = "in",
    limitsize = FALSE
  )
}
end <- Sys.time()

log_debug("Script finished in", format(end - start))
