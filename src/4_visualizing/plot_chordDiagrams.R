source("r/log_debug.R")
log_debug("This script plots the interactive chord diagrams (for the SI)")

start <- Sys.time()
library(chorddiag)
library(data.table)
library(dplyr)
library(plotly)
library(splitstackshape)
library(tibble)
library(tidyr)
source(file = "paths.R")
source("r/draw_chord.R")
source("r/palettes.R")
source("r/y_as_na.R")

pairs_metadata <-
  read_delim(
    file = file.path(
      pathDataProcessed,
      pathLastFrozen
    ),
    col_types = cols(.default = "c"),
    locale = locales
  ) %>%
  data.table() %>%
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
  mutate_all(as.character)

pairs_metadata[] <-
  lapply(pairs_metadata, function(x) {
    y_as_na(x, "")
  })

log_debug("... drawing big chord diagram")
top_big_chord_bio <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_01domain) &
      !is.na(structure_taxonomy_npclassifier_01pathway)
  ) %>%
  group_by(organism_taxonomy_01domain) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_01domain) %>%
  head(6)

top_big_chord_bio <- top_big_chord_bio$organism_taxonomy_01domain

top_big_chord_chemo <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_01domain) &
      !is.na(structure_taxonomy_npclassifier_01pathway)
  ) %>%
  filter(organism_taxonomy_01domain %in% top_big_chord_bio) %>%
  group_by(structure_taxonomy_npclassifier_01pathway) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_01pathway) %>%
  head(12)

top_big_chord_chemo <-
  top_big_chord_chemo$structure_taxonomy_npclassifier_01pathway

chord_big <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_01domain",
  chemical_level = "structure_taxonomy_npclassifier_01pathway",
  chemical_filter_value = top_big_chord_chemo,
  chemical_filter_level = "structure_taxonomy_npclassifier_01pathway",
  biological_filter_value = top_big_chord_bio,
  biological_filter_level = "organism_taxonomy_01domain",
  palette = paired_palette_med
)

log_debug("... drawing medium chord diagram")
top_organism_med <- pairs_metadata %>%
  filter(!is.na(organism_taxonomy_06family)) %>%
  filter(structure_taxonomy_npclassifier_01pathway == "Alkaloids") %>%
  group_by(organism_taxonomy_06family) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_06family) %>%
  head(12)

top_organism_med <- top_organism_med$organism_taxonomy_06family

top_chemo_med <- pairs_metadata %>%
  filter(!is.na(organism_taxonomy_06family)) %>%
  filter(structure_taxonomy_npclassifier_01pathway == "Alkaloids") %>%
  filter(!is.na(structure_taxonomy_npclassifier_02superclass)) %>%
  filter(organism_taxonomy_06family %in% top_organism_med) %>%
  group_by(structure_taxonomy_npclassifier_02superclass) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_02superclass) %>%
  head(18)

top_chemo_med <-
  top_chemo_med$structure_taxonomy_npclassifier_02superclass

chord_med <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_06family",
  chemical_level = "structure_taxonomy_npclassifier_02superclass",
  chemical_filter_value = top_chemo_med,
  chemical_filter_level = "structure_taxonomy_npclassifier_02superclass",
  biological_filter_value = top_organism_med,
  biological_filter_level = "organism_taxonomy_06family",
  palette = paired_palette_30
)

log_debug("... drawing small chord diagram")
top_organism_sma <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_08genus == "Erythroxylum") %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(12)

top_organism_sma <- top_organism_sma$organism_taxonomy_09species

top_chemo_sma <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_08genus == "Erythroxylum") %>%
  filter(organism_taxonomy_09species %in% top_organism_sma) %>%
  filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class) %>%
  head(12)

top_chemo_sma <-
  top_chemo_sma$structure_taxonomy_npclassifier_03class

chord_sma <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_sma,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_sma,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_big
)

log_debug("... drawing Ranunculaceae chord diagram")
top_organism_ranunculaceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Ranunculaceae") %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(12)

top_organism_ranunculaceae <-
  top_organism_ranunculaceae$organism_taxonomy_09species

top_chemo_ranunculaceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Ranunculaceae") %>%
  filter(organism_taxonomy_09species %in% top_organism_ranunculaceae) %>%
  filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class) %>%
  head(12)

top_chemo_ranunculaceae <-
  top_chemo_ranunculaceae$structure_taxonomy_npclassifier_03class

chord_ranunculaceae <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_ranunculaceae,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_ranunculaceae,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_big
)

log_debug("... drawing Papaveraceae chord diagram")
top_organism_papaveraceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Papaveraceae") %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(12)

top_organism_papaveraceae <-
  top_organism_papaveraceae$organism_taxonomy_09species

top_chemo_papaveraceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Papaveraceae") %>%
  filter(organism_taxonomy_09species %in% top_organism_papaveraceae) %>%
  filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class) %>%
  head(12)

top_chemo_papaveraceae <-
  top_chemo_papaveraceae$structure_taxonomy_npclassifier_03class

chord_papaveraceae <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_papaveraceae,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_papaveraceae,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_big
)

log_debug("... drawing Gentianaceae chord diagram")
top_organism_gentianaceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Gentianaceae") %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(36)

top_organism_gentianaceae <-
  top_organism_gentianaceae$organism_taxonomy_09species

top_chemo_gentianaceae <- pairs_metadata %>%
  filter(
    !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  filter(organism_taxonomy_06family == "Gentianaceae") %>%
  filter(structure_taxonomy_npclassifier_02superclass == "Xanthones") %>%
  filter(organism_taxonomy_09species %in% top_organism_gentianaceae) %>%
  filter(!is.na(structure_inchikey)) %>%
  group_by(structure_inchikey) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_inchikey) %>%
  head(144)

top_chemo_gentianaceae <-
  top_chemo_gentianaceae$structure_inchikey

chord_gentianaceae <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_inchikey",
  chemical_filter_value = top_chemo_gentianaceae,
  chemical_filter_level = "structure_inchikey",
  biological_filter_value = top_organism_gentianaceae,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_big
)

log_debug("... drawing top N chord diagrams ...")
log_debug("... top 06")
top_organism_06 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(6)

top_organism_06 <-
  top_organism_06$organism_taxonomy_09species

top_chemo_06 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(structure_taxonomy_npclassifier_03class)
  ) %>%
  filter(organism_name %in% top_organism_06) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class, .keep_all = TRUE) %>%
  head(6)

top_chemo_06 <-
  top_chemo_06$structure_taxonomy_npclassifier_03class

chord_06 <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_06,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_06,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_sma
)

log_debug("... top 12")
top_organism_12 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(12)

top_organism_12 <-
  top_organism_12$organism_taxonomy_09species

top_chemo_12 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(structure_taxonomy_npclassifier_03class)
  ) %>%
  filter(organism_name %in% top_organism_12) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class, .keep_all = TRUE) %>%
  head(12)

top_chemo_12 <-
  top_chemo_12$structure_taxonomy_npclassifier_03class

chord_12 <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_12,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_12,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_big
)

log_debug("... top 24")
top_organism_24 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(organism_taxonomy_09species)
  ) %>%
  group_by(organism_taxonomy_09species) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(organism_taxonomy_09species) %>%
  head(24)

top_organism_24 <-
  top_organism_24$organism_taxonomy_09species

top_chemo_24 <- pairs_metadata %>%
  filter(
    !is.na(organism_taxonomy_09species) &
      !is.na(structure_taxonomy_npclassifier_03class) &
      !is.na(structure_taxonomy_npclassifier_03class)
  ) %>%
  filter(organism_name %in% top_organism_24) %>%
  group_by(structure_taxonomy_npclassifier_03class) %>%
  add_count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  distinct(structure_taxonomy_npclassifier_03class, .keep_all = TRUE) %>%
  head(24)

top_chemo_24 <-
  top_chemo_24$structure_taxonomy_npclassifier_03class

chord_24 <- draw_chord(
  data = pairs_metadata,
  biological_level = "organism_taxonomy_09species",
  chemical_level = "structure_taxonomy_npclassifier_03class",
  chemical_filter_value = top_chemo_24,
  chemical_filter_level = "structure_taxonomy_npclassifier_03class",
  biological_filter_value = top_organism_24,
  biological_filter_level = "organism_taxonomy_09species",
  palette = paired_palette_meg
)

if (mode == "full") {
  create_dir("../res/html")
  setwd("../res/html")
  htmlwidgets::saveWidget(
    widget = as_widget(chord_big),
    file = "chord_big.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_med),
    file = "chord_med.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_sma),
    file = "chord_sma.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_ranunculaceae),
    file = "chord_ranunculaceae.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_papaveraceae),
    file = "chord_papaveraceae.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_gentianaceae),
    file = "chord_gentianaceae.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_06),
    file = "chord_06.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_12),
    file = "chord_12.html"
  )

  htmlwidgets::saveWidget(
    widget = as_widget(chord_24),
    file = "chord_24.html"
  )
}
end <- Sys.time()

log_debug("Script finished in", format(end - start))
