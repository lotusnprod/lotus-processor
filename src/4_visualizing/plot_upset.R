source("r/log_debug.R")
log_debug("This script plots the upset plots ...")

start <- Sys.time()

log_debug("sourcing ...")
log_debug("... paths")
source("paths.R")

if (mode == "full") {
  log_debug("... libraries")
  library(data.table)
  library(dplyr)
  library(splitstackshape)
  library(readr)
  library(tidyr)
  library(UpSetR)
  source("r/getGraphChemicalClass.R")
  source("r/getGraphStudiedPlant.R")
  source("r/prepare_upset.R")
  source("r/y_as_na.R")

  log_debug("loading files, this may take a while ...")

  log_debug("... pretty names")
  prettyNames <-
    read_delim(file = "../docs/prettyDBNames.tsv")

  log_debug("... open DB")
  openDb <-
    read_delim(
      file = pathDataInterimTablesAnalyzedPlatinum,
      col_types = cols(.default = "c")
    ) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey,
      referenceCleanedTitle,
      .keep_all = TRUE
    ) %>%
    filter(database != "wikidata")

  log_debug("... closed DBs")
  closedDb <-
    read_delim(file = file.path(pathDataInterimTablesAnalyzed, "closed.tsv.gz"),
               col_types = cols(.default = "c")) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey,
      structureCleaned_inchikey2D,
      .keep_all = TRUE
    )

  log_debug("... metadata ...")
  pairs_metadata <- fread(
    file = file.path(
      pathDataProcessed,
      pathLastFrozen
    )
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

  log_debug("joining closed DBs and open DB")
  inhouseDb <- bind_rows(closedDb, openDb) %>%
    left_join(., prettyNames) %>%
    select(-database) %>%
    select(
      database = prettyDataBase,
      everything()
    )

  log_debug("counting sub-DBs")
  dbNum <-
    as.numeric(nrow(inhouseDb %>%
      distinct(inhouseDb$database)))

  log_debug("drawing upset plot of structures repartition in sub-DBs")
  inhouseDb_structures_2plot_wide <-
    prepare_upset(
      table = inhouseDb,
      group = "database",
      variable = "structureCleaned_inchikey2D"
    )
  pdf(
    file = file.path("../res", "upset_structures.pdf"),
    width = 16,
    height = 9
  )
  upset(
    inhouseDb_structures_2plot_wide,
    nsets = 10,
    order.by = "freq",
    nintersects = 20,
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 1.6,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per database",
    set_size.scale_max = 250000,
    set_size.show = TRUE
  )
  dev.off()

  log_debug("drawing upset plot of organisms repartition in sub-DBs")
  inhouseDb_organism_2plot_wide <-
    prepare_upset(
      table = inhouseDb,
      group = "database",
      variable = "organismCleaned"
    )
  pdf(
    file = file.path("../res", "upset_organisms.pdf"),
    width = 16,
    height = 9
  )
  upset(
    inhouseDb_organism_2plot_wide,
    nsets = 10,
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 20,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique organisms per intersection",
    sets.x.label = "Unique organisms per database",
    set_size.scale_max = 45000,
    set_size.show = TRUE
  )
  dev.off()

  log_debug("drawing upset plot of pairs repartition in sub-DBs")
  inhouseDbPairs_2plot_wide <-
    prepare_upset(
      table = inhouseDb,
      group = "database",
      variable = c("structureCleaned_inchikey2D", "organismCleaned")
    )
  pdf(
    file = file.path("../res", "upset_pairs.pdf"),
    width = 16,
    height = 9
  )
  upset(
    inhouseDbPairs_2plot_wide,
    nsets = 10,
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 20,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 1.6,
    mainbar.y.label = "Unique pairs per intersection",
    sets.x.label = "Unique pairs per database",
    set_size.show = TRUE,
    set_size.scale_max = 300000
  )
  dev.off()

  log_debug("adding metadata for more detailed analyzis ...")
  log_debug("... classyfire")
  inhouseDbMeta <- left_join(
    inhouseDb %>% distinct(
      database,
      organism_name = organismCleaned,
      structure_inchikey = structureCleanedInchikey,
      structure_inchikey_2D = structureCleaned_inchikey2D
    ),
    pairs_metadata %>% distinct(organism_name, structure_inchikey, .keep_all = TRUE)
  )

  chemo <- inhouseDbMeta %>%
    filter(!is.na(structure_taxonomy_npclassifier_01pathway)) %>%
    distinct(
      structure_inchikey_2D,
      structure_taxonomy_npclassifier_01pathway
    )

  chemo3D <- inhouseDbMeta %>%
    filter(!is.na(structure_taxonomy_npclassifier_01pathway)) %>%
    distinct(
      structure_inchikey,
      structure_taxonomy_npclassifier_01pathway
    )

  bio <- inhouseDbMeta %>%
    filter(!is.na(organism_taxonomy_03phylum)) %>%
    distinct(
      organism_name,
      organism_taxonomy_03phylum
    )

  log_debug("drawing upset plot of stigmastenol repartition in sub-DBs")
  inhouseDb_most_structures <- inhouseDbMeta %>%
    filter(!is.na(structure_inchikey_2D)) %>%
    count(structure_inchikey_2D) %>%
    arrange(desc(n))

  mostinchi <- as.character(inhouseDb_most_structures[1, 1])

  inhouseDb_most_structures_2plot <- inhouseDbMeta %>%
    filter(structure_inchikey_2D == mostinchi) %>%
    filter(!is.na(organism_taxonomy_03phylum))

  inhouseDb_most_structures_2plot_wide <-
    prepare_upset_complex(
      table = inhouseDb_most_structures_2plot,
      group = c("database", "organism_name"),
      ## order is important
      variable = "structure_inchikey_2D"
    ) %>%
    left_join(
      .,
      bio
    ) %>%
    distinct(structure_inchikey_2D,
      organism_name,
      .keep_all = TRUE
    )

  mostkingdom <- inhouseDb_most_structures_2plot_wide %>%
    filter(!is.na(organism_taxonomy_03phylum)) %>%
    count(organism_taxonomy_03phylum) %>%
    arrange(desc(n))

  dbnumostinchi <- as.numeric(nrow(
    inhouseDbMeta %>%
      filter(structure_inchikey_2D == mostinchi) %>%
      distinct(database)
  ))

  pdf(
    file = file.path("../res", "upset_stigmastenol.pdf"),
    width = 16,
    height = 9
  )

  upset(
    inhouseDb_most_structures_2plot_wide,
    nsets = 10,
    query.legend = "top",
    queries = list(
      list(
        query = elements,
        params = list(
          "organism_taxonomy_03phylum",
          c(
            mostkingdom[1, 1],
            mostkingdom[2, 1],
            mostkingdom[3, 1]
          )
        ),
        active = TRUE,
        color = "#b2df8a",
        query.name = mostkingdom[3, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_03phylum",
          c(
            mostkingdom[1, 1],
            mostkingdom[2, 1]
          )
        ),
        active = TRUE,
        color = "#1f78b4",
        query.name = mostkingdom[2, 1]
      ),
      list(
        query = elements,
        params = list(
          "organism_taxonomy_03phylum",
          mostkingdom[1, 1]
        ),
        active = TRUE,
        color = "#a6cee3",
        query.name = mostkingdom[1, 1]
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
    mainbar.y.label = "Unique organisms per intersection",
    sets.x.label = "Unique organisms per database",
    set_size.show = TRUE,
    set_size.scale_max = 2500
  )
  dev.off()

  log_debug("getting list of most studied plants for exploration")
  inhouseDb_most_plant <- inhouseDbMeta %>%
    filter(!is.na(organism_name)) %>%
    filter(organism_taxonomy_02kingdom == "Archaeplastida") %>%
    count(organism_name) %>%
    arrange(desc(n))

  log_debug("... drawing Cannabis sativa metabolites repartition")
  pdf(
    file = file.path("../res", "upset_cannabisSativa.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Cannabis sativa")
  dev.off()

  log_debug("... drawing Tripterygium wilfordii metabolites repartition")
  pdf(
    file = file.path("../res", "upset_tripterygiumWilfordii.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Tripterygium wilfordii")
  dev.off()

  log_debug("... drawing Citrus reticulata metabolites repartition")
  pdf(
    file = file.path("../res", "upset_citrusReticulata.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Citrus reticulata")
  dev.off()

  log_debug("... drawing Camellia sinensis metabolites repartition")
  pdf(
    file = file.path("../res", "upset_camelliaSinensis.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Camellia sinensis")
  dev.off()

  log_debug("... drawing Panax ginseng metabolites repartition")
  pdf(
    file = file.path("../res", "upset_panaxGinseng.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Panax ginseng")
  dev.off()

  log_debug("... drawing Vitis vinifera metabolites repartition")
  pdf(
    file = file.path("../res", "upset_vitisVinifera.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Vitis vinifera")
  dev.off()

  log_debug("... drawing Zingiber officinale metabolites repartition")
  pdf(
    file = file.path("../res", "upset_zingiberOfficinale.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Zingiber officinale")
  dev.off()

  log_debug("... drawing Capsicum annuum metabolites repartition")
  pdf(
    file = file.path("../res", "upset_capsicumAnnuum.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Capsicum annuum")
  dev.off()

  log_debug("... drawing Arabidopsis thaliana metabolites repartition")
  pdf(
    file = file.path("../res", "upset_arabidopsisThaliana.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Arabidopsis thaliana")
  dev.off()

  log_debug("getting list of most studied chemical subclasses for exploration")
  inhouseDb_most_chemical_subclass <- inhouseDbMeta %>%
    filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
    count(structure_taxonomy_npclassifier_03class) %>%
    arrange(desc(n))

  log_debug("... drawing flavonols repartition")
  pdf(
    file = file.path("../res", "upset_flavonols.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Flavonols")
  dev.off()

  log_debug("... drawing Iridoids monoterpenoids repartition")
  pdf(
    file = file.path("../res", "upset_iridoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Iridoids monoterpenoids")
  dev.off()

  log_debug("... drawing Limonoids repartition")
  pdf(
    file = file.path("../res", "upset_limonoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Limonoids")
  dev.off()

  log_debug("... drawing Quassinoids repartition")
  pdf(
    file = file.path("../res", "upset_quassinoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Quassinoids")
  dev.off()
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
