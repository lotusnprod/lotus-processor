cat("This script plots the upset plots ... \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

if (mode == "full") {
  cat("... libraries \n")
  library(data.table)
  library(splitstackshape)
  library(tidyverse)
  library(UpSetR)
  source("r/vroom_safe.R")
  source("r/getGraphChemicalClass.R")
  source("r/getGraphStudiedPlant.R")
  source("r/prepare_upset.R")
  source("r/y_as_na.R")

  cat("loading files, this may take a while ... \n")

  cat("... pretty names \n")
  prettyNames <-
    vroom_read_safe(path = "../docs/prettyDBNames.tsv")

  cat("... open DB \n")
  openDb <-
    vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey,
      referenceCleanedTitle,
      .keep_all = TRUE
    ) %>%
    tibble()

  cat("... open DB \n")
  openDbMaximal <-
    vroom_read_safe(path = pathDataInterimTablesCuratedTableMaximal) %>%
    tibble()

  cat("... DNP DB \n")
  dnpDb <-
    vroom_read_safe(path = file.path(pathDataInterimTablesAnalysed, "dnp.tsv.gz")) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey,
      structureCleaned_inchikey2D,
      .keep_all = TRUE
    ) %>%
    tibble()

  cat("... metadata ... \n")
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
    mutate_all(as.character) %>%
    tibble()

  pairs_metadata[] <-
    lapply(pairs_metadata, function(x) {
      y_as_na(x, "")
    })

  cat("joining DNP and openDB \n")
  inhouseDb <- bind_rows(dnpDb, openDb) %>%
    left_join(., prettyNames) %>%
    select(-database) %>%
    select(
      database = prettyDataBase,
      everything()
    )

  cat("counting sub-DBs \n")
  dbNum <-
    as.numeric(nrow(inhouseDb %>%
      distinct(inhouseDb$database)))

  cat("drawing upset plot of structures repartition in sub-DBs \n")
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

  cat("drawing upset plot of organisms repartition in sub-DBs \n")
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

  cat("drawing upset plot of pairs repartition in sub-DBs \n")
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

  cat("adding metadata for more detailed analysis ... \n")
  cat("... classyfire \n")
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

  cat("drawing upset plot of stigmastenol repartition in sub-DBs \n")
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

  cat("getting list of most studied plants for exploration \n")
  inhouseDb_most_plant <- inhouseDbMeta %>%
    filter(!is.na(organism_name)) %>%
    filter(organism_taxonomy_02kingdom == "Archaeplastida") %>%
    count(organism_name) %>%
    arrange(desc(n))

  cat("... drawing Cannabis sativa metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_cannabisSativa.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Cannabis sativa")
  dev.off()

  cat("... drawing Tripterygium wilfordii metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_tripterygiumWilfordii.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Tripterygium wilfordii")
  dev.off()

  cat("... drawing Citrus aurantium metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_citrusAurantium.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Citrus aurantium")
  dev.off()

  cat("... drawing Camellia sinensis metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_camelliaSinensis.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Camellia sinensis")
  dev.off()

  cat("... drawing Panax ginseng metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_panaxGinseng.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Panax ginseng")
  dev.off()

  cat("... drawing Vitis vinifera metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_vitisVinifera.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Vitis vinifera")
  dev.off()

  cat("... drawing Zingiber officinale metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_zingiberOfficinale.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Zingiber officinale")
  dev.off()

  cat("... drawing Capsicum annuum metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_capsicumAnnuum.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Capsicum annuum")
  dev.off()

  cat("... drawing Arabidopsis thaliana metabolites repartition \n")
  pdf(
    file = file.path("../res", "upset_arabidopsisThaliana.pdf"),
    width = 16,
    height = 9
  )
  getGraphStudiedPlant(plant = "Arabidopsis thaliana")
  dev.off()

  cat("getting list of most studied chemical subclasses for exploration \n")
  inhouseDb_most_chemical_subclass <- inhouseDbMeta %>%
    filter(!is.na(structure_taxonomy_npclassifier_03class)) %>%
    count(structure_taxonomy_npclassifier_03class) %>%
    arrange(desc(n))

  cat("... drawing flavonols repartition \n")
  pdf(
    file = file.path("../res", "upset_flavonols.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Flavonols")
  dev.off()

  cat("... drawing Iridoids monoterpenoids repartition \n")
  pdf(
    file = file.path("../res", "upset_iridoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Iridoids monoterpenoids")
  dev.off()

  cat("... drawing Limonoids repartition \n")
  pdf(
    file = file.path("../res", "upset_limonoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Limonoids")
  dev.off()

  cat("... drawing Quassinoids repartition \n")
  pdf(
    file = file.path("../res", "upset_quassinoids.pdf"),
    width = 16,
    height = 9
  )
  getGraphChemicalClass(subclass = "Quassinoids")
  dev.off()
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")