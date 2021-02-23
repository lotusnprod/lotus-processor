cat("This script plots the wonderful tree. \n")
cat("It currently needs 'tree_prepareTree.R' to be run before. \n")

start <- Sys.time()

library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
source(file = "r/colors.R")

test <- pairs_metadata %>%
  filter(Family %in% families_restricted$Family) %>%
  # filter(!is.na(`03class`)) %>%
  filter(!is.na(class)) %>%
  distinct(Family,
    structureCleaned_inchikey2D,
    .keep_all = TRUE
  ) %>%
  group_by(Family) %>%
  add_count(name = "n") %>%
  # group_by(`03class`) %>%
  group_by(class) %>%
  add_count(name = "m") %>%
  group_by(
    Family,
    Genus
  ) %>%
  add_count(name = "o") %>%
  # group_by(Family, `03class`) %>%
  group_by(Family, class) %>%
  add_count(name = "p") %>%
  mutate(q = p^2 / (m * n)) %>%
  # distinct(Family,
  #   `03class`,
  #   q,
  #   o,
  #   .keep_all = TRUE
  # ) %>%
  distinct(Family,
    class,
    q,
    o,
    .keep_all = TRUE
  ) %>%
  group_by(Family) %>%
  filter(q == max(q)) %>%
  ungroup() %>%
  distinct(Family, q, .keep_all = TRUE)

tr_restricted$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_restricted$tip.label
  )

taxonomy <-
  left_join(
    families_restricted,
    pairs_metadata %>%
      distinct(
        Domain,
        Kingdom,
        Phylum,
        Class,
        Order,
        Infraorder,
        Family
      )
  )

info <- taxonomy %>%
  select(
    id = Family,
    everything()
  ) %>%
  mutate(Kingdom = fct_reorder(Kingdom, !is.na(Domain))) %>%
  mutate(Phylum = fct_reorder(Phylum, !is.na(Kingdom)))

bar_data <- test %>%
  select(
    id = Family,
    everything(),
    -Domain,
    -Kingdom,
    -Phylum,
    -Class,
    -Order
  ) %>%
  # mutate(`02superclass` = fct_reorder(`02superclass`, `01kingdom`)) %>%
  # mutate(`03class` = fct_reorder(`03class`, !is.na(`02superclass`)))
  mutate(superclass = fct_reorder(superclass, pathway)) %>%
  mutate(class = fct_reorder(class, !is.na(superclass)))

bar_data_temp <- bar_data

info_temp <- info %>%
  filter(id %in% bar_data_temp$id)

ott_temp <-
  ott_in_tree_restricted[names(ott_in_tree_restricted) %in% bar_data_temp$id]

ott_temp <- ott_temp[!duplicated(ott_temp)]

tr_temp <- tol_induced_subtree(ott_ids = ott_temp)

tr_temp$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_temp$tip.label
  )

p <- ggtree(tr = tr_temp, layout = "circular")

p <- p %<+%
  info_temp +
  # geom_label(aes(x=branch, label=Order)) +
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
        # fill = `02superclass`
        fill = pathway
      ),
    offset = rel(0.2),
    pwidth = rel(1.6),
    orientation = "y",
    stat = "identity",
  ) +
  geom_tiplab(
    mapping = aes(
      # label = paste(`02superclass`, `03class`, sep = " - "),
      label = paste(superclass, class, sep = " - "),
      x = 26 + 63 * q ## to adapt with the data
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
  # scale_fill_manual(
  # values = strsplit(
  # x = paired[4:18],
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

p

ggsave(
  filename = file.path(pathDataProcessedFigures, "magicTree.pdf"),
  plot = p,
  width = 100,
  height = 100,
  units = "in",
  limitsize = FALSE
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
