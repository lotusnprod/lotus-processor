cat("This script plots the interactive chord diagrams (for the SI) \n")

start <- Sys.time()
library(chorddiag)
library(data.table)
library(plotly)
library(splitstackshape)
library(tidyverse)
source(file = "paths.R")
source("r/draw_chord.R")
source("r/palettes.R")
source("r/y_as_na.R")

pairs_metadata <- fread(
  input = file.path(
    pathDataProcessed,
    "210223_frozen_metadata.csv.gz"
  ),
  na.strings = ""
) %>%
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

cat("... drawing big chord diagram \n")
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

setwd("../res/html")
htmlwidgets::saveWidget(
  widget = as_widget(chord_big),
  file = "chord_big.html"
)

cat("... drawing medium chord diagram \n")
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_med),
  file = "chord_med.html"
)

cat("... drawing small chord diagram \n")
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_sma),
  file = "chord_sma.html"
)

cat("... drawing Ranunculaceae chord diagram \n")
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_ranunculaceae),
  file = "chord_ranunculaceae.html"
)

cat("... drawing Papaveraceae chord diagram \n")
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_papaveraceae),
  file = "chord_papaveraceae.html"
)

cat("... drawing Gentianaceae chord diagram \n")
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_gentianaceae),
  file = "chord_gentianaceae.html"
)

cat("... drawing top N chord diagrams ... \n")
cat("... top 06 \n")
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
  distinct(structure_taxonomy_npclassifier_03class,
    .keep_all = TRUE
  ) %>%
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_06),
  file = "chord_06.html"
)

cat("... top 12 \n")
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
  distinct(structure_taxonomy_npclassifier_03class,
    .keep_all = TRUE
  ) %>%
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_12),
  file = "chord_12.html"
)

cat("... top 24 \n")
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
  distinct(structure_taxonomy_npclassifier_03class,
    .keep_all = TRUE
  ) %>%
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

htmlwidgets::saveWidget(
  widget = as_widget(chord_24),
  file = "chord_24.html"
)

cat("Also exporting static \n")

orca(
  p = chord_24,
  file = file.path("..", "chord_24.pdf")
)

orca(
  p = chord_12,
  file = file.path("..", "chord_12.pdf")
)

orca(
  p = chord_06,
  file = file.path("..", "chord_06.pdf")
)

orca(
  p = chord_big,
  file = file.path("..", "chord_big.pdf")
)

orca(
  p = chord_med,
  file = file.path("..", "chord_med.pdf")
)

orca(
  p = chord_sma,
  file = file.path("..", "chord_sma.pdf")
)

orca(
  p = chord_ranunculaceae,
  file = file.path("..", "chord_ranunculaceae.pdf")
)

orca(
  p = chord_papaveraceae,
  file = file.path("..", "chord_papaveraceae.pdf")
)

orca(
  p = chord_gentianaceae,
  file = file.path("..", "chord_gentianaceae.pdf")
)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
