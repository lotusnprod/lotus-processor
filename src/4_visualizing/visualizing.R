cat("This script draws figures associated with DBs content \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")
cat("... functions \n")
source("functions.R")

if (mode != "test") {
  cat("loading files, if running fullmode, this may take a while ... \n")
  cat("... open DB \n")
  openDb <- read_delim(
    file = gzfile(pathDataInterimTablesAnalysedPlatinum),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey3D,
      referenceCleanedTitle,
      .keep_all = TRUE
    ) %>%
    tibble()

  cat("... open DB \n")
  openDbMaximal <- read_delim(
    file = gzfile(pathDataInterimTablesCuratedTableMaximal),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble()

  cat("... DNP DB \n")
  dnpDb <- read_delim(
    file = gzfile(file.path(
      pathDataInterimTablesAnalysed, "dnp.tsv.gz"
    )),
    col_types = cols(.default = "c"),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey3D,
      structureCleanedInchikey2D,
      .keep_all = TRUE
    ) %>%
    tibble()

  cat("... metadata ... \n")
  cat("... organisms \n")
  organismMetadata <- read_delim(
    file = gzfile(pathDataInterimDictionariesOrganismMetadata),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble() %>%
    mutate(
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Viridiplantae",
        yes = "Plantae",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Metazoa",
        yes = "Animalia",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Proteobacteria",
        yes = "Bacteria",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Protozoa",
        yes = "Protista",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_2phylum = ifelse(
        test = organismCleaned_dbTaxo_2phylum == "Streptophyta" |
          organismCleaned_dbTaxo_2phylum == "Magnoliophyta" |
          organismCleaned_dbTaxo_2phylum == "Gymnospermophyta",
        yes = "Tracheophyta",
        no = organismCleaned_dbTaxo_2phylum
      ),
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Hepaticae",
        yes = "Plantae",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_2phylum = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Hepaticae",
        yes = "Marchantiophyta",
        no = organismCleaned_dbTaxo_2phylum
      ),
      organismCleaned_dbTaxo_1kingdom = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Cyanobacteria",
        yes = "Bacteria",
        no = organismCleaned_dbTaxo_1kingdom
      ),
      organismCleaned_dbTaxo_2phylum = ifelse(
        test = organismCleaned_dbTaxo_1kingdom == "Cyanobacteria",
        yes = "Cyanobacteria",
        no = organismCleaned_dbTaxo_2phylum
      )
    )

  cat("... structures metadata \n")
  structureMetadata_1 <- read_delim(
    file = gzfile(pathDataInterimDictionariesStructureMetadata),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    distinct(
      structureCleanedSmiles,
      structureCleanedInchi,
      structureCleanedInchikey3D,
      structureCleaned_inchikey2D,
      structureCleaned_stereocenters_unspecified
    ) %>%
    tibble()

  cat("... structures classification \n")
  structureMetadata_2 <- read_delim(
    file = gzfile(
      "../data/interim/09_structureContextualizer/classy.tsv.gz"
    ),
    # dirty residue, will have to change
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
    tibble()

  cat("joining DNP and openDB \n")
  inhouseDb <-
    bind_rows(dnpDb, openDb)

  cat("counting sub-DBs \n")
  dbNum <-
    as.numeric(nrow(inhouseDb %>%
      distinct(inhouseDb$database)))

  cat("ensuring directories exist ... \n")
  ifelse(
    test = !dir.exists(pathDataProcessed),
    yes = dir.create(pathDataProcessed),
    no = paste(pathDataProcessed, "exists")
  )

  ifelse(
    test = !dir.exists(pathDataProcessedFigures),
    yes = dir.create(pathDataProcessedFigures),
    no = paste(pathDataProcessedFigures, "exists")
  )

  ifelse(
    test = !dir.exists(pathDataProcessedFiguresHtml),
    yes = dir.create(pathDataProcessedFiguresHtml),
    no = paste(pathDataProcessedFiguresHtml, "exists")
  )

  cat("drawing upset plot of structures repartition in sub-DBs \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "structures.pdf"),
    width = 16,
    height = 9
  )

  inhouseDb_structures_2plot <- inhouseDb %>%
    distinct(structureCleanedInchikey2D, database) %>%
    group_by(database) %>%
    count(structureCleanedInchikey2D) %>%
    ungroup()

  inhouseDb_structures_2plot_wide <- inhouseDb_structures_2plot %>%
    pivot_wider(
      names_from = database,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    distinct(structureCleanedInchikey2D, .keep_all = TRUE) %>%
    data.frame()

  upset(
    inhouseDb_structures_2plot_wide,
    nsets = 10,
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 20,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per database",
    set_size.show = TRUE
  )
  dev.off()

  cat("drawing upset plot of organisms repartition in sub-DBs \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "organisms.pdf"),
    width = 16,
    height = 9
  )

  inhouseDb_organism_2plot <- inhouseDb %>%
    filter(!is.na(organismCleaned)) %>%
    distinct(organismCleaned, database) %>%
    group_by(database) %>%
    count(organismCleaned) %>%
    ungroup()

  inhouseDb_organism_2plot_wide <-
    inhouseDb_organism_2plot %>%
    pivot_wider(
      names_from = database,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    distinct(organismCleaned, .keep_all = TRUE) %>%
    data.frame()


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
    mainbar.y.label = "Unique biological sources per intersection",
    sets.x.label = "Unique sources per database",
    set_size.show = TRUE
  )
  dev.off()

  cat("drawing upset plot of pairs repartition in sub-DBs \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "pairs.pdf"),
    width = 16,
    height = 9
  )

  inhouseDbPairs_2plot <- inhouseDb %>%
    filter(!is.na(structureCleanedInchikey2D) &
      !is.na(organismCleaned)) %>%
    distinct(structureCleanedInchikey2D,
      organismCleaned,
      database,
      .keep_all = TRUE
    ) %>%
    group_by(database) %>%
    count(structureCleanedInchikey2D, organismCleaned) %>%
    ungroup()

  inhouseDbPairs_2plot_wide <- inhouseDbPairs_2plot %>%
    pivot_wider(
      names_from = database,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(3:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(3:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    distinct(structureCleanedInchikey2D, organismCleaned, .keep_all = TRUE) %>%
    data.frame()

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
    text.scale = 2,
    mainbar.y.label = "Unique pairs per intersection",
    sets.x.label = "Unique pairs per database",
    set_size.show = TRUE
  )
  dev.off()

  cat("adding metadata for more detailed analysis ... \n")
  cat("... inhouse DB \n")
  inhouseDbMeta <- left_join(inhouseDb, structureMetadata_2)

  inhouseDbMeta <- left_join(inhouseDbMeta, organismMetadata) %>%
    arrange(desc(organismCleaned_dbTaxo_8variety)) %>%
    arrange(desc(organismCleaned_dbTaxo_7species)) %>%
    arrange(desc(organismCleaned_dbTaxo_6genus)) %>%
    arrange(desc(organismCleaned_dbTaxo_5family)) %>%
    arrange(desc(organismCleaned_dbTaxo_4order)) %>%
    arrange(desc(organismCleaned_dbTaxo_3class)) %>%
    arrange(desc(organismCleaned_dbTaxo_2phylum)) %>%
    arrange(desc(organismCleaned_dbTaxo_1kingdom)) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey2D,
      referenceCleanedTitle,
      .keep_all = TRUE
    )

  cat("... open DB \n")
  openDbMeta <- left_join(openDb, structureMetadata_2)

  openDbMeta <- left_join(openDbMeta, organismMetadata) %>%
    arrange(desc(organismCleaned_dbTaxo_8variety)) %>%
    arrange(desc(organismCleaned_dbTaxo_7species)) %>%
    arrange(desc(organismCleaned_dbTaxo_6genus)) %>%
    arrange(desc(organismCleaned_dbTaxo_5family)) %>%
    arrange(desc(organismCleaned_dbTaxo_4order)) %>%
    arrange(desc(organismCleaned_dbTaxo_3class)) %>%
    arrange(desc(organismCleaned_dbTaxo_2phylum)) %>%
    arrange(desc(organismCleaned_dbTaxo_1kingdom)) %>%
    distinct(
      database,
      organismCleaned,
      structureCleanedInchikey2D,
      referenceCleanedTitle,
      .keep_all = TRUE
    )

  chemo <- inhouseDbMeta %>%
    filter(!is.na(structureCleaned_2superclass)) %>%
    distinct(
      structureCleanedInchikey2D,
      structureCleaned_2superclass
    )

  chemo3D <- inhouseDbMeta %>%
    filter(!is.na(structureCleaned_2superclass)) %>%
    distinct(
      structureCleanedInchikey3D,
      structureCleaned_2superclass
    )

  bio <- inhouseDbMeta %>%
    filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
    distinct(
      organismCleaned,
      organismCleaned_dbTaxo_1kingdom
    )

  try({
    cat("drawing upset plot of stigmastenol repartition in sub-DBs \n")
    pdf(
      file = file.path(pathDataProcessedFigures, "stigmastenol.pdf"),
      width = 16,
      height = 9
    )
    inhouseDb_most_structures <- inhouseDbMeta %>%
      filter(!is.na(structureCleanedInchikey2D)) %>%
      count(structureCleanedInchikey2D) %>%
      arrange(desc(n))
    mostinchi <- as.character(inhouseDb_most_structures[1, 1])
    inhouseDb_most_structures_2plot <-
      inhouseDbMeta %>%
      filter(structureCleanedInchikey2D == mostinchi) %>%
      filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        database,
        .keep_all = TRUE
      ) %>%
      group_by(database, organismCleaned) %>%
      count(structureCleanedInchikey2D) %>%
      ungroup()
    inhouseDb_most_structures_2plot_wide <-
      inhouseDb_most_structures_2plot %>%
      pivot_wider(
        names_from = database,
        values_from = n
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = is.na(.),
          values = 0
        )
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = . >= 1,
          values = 1
        )
      ) %>%
      distinct(organismCleaned, .keep_all = TRUE) %>%
      data.frame()
    inhouseDb_most_structures_2plot_wide <-
      left_join(
        inhouseDb_most_structures_2plot_wide,
        bio
      ) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        .keep_all = TRUE
      )
    mostkingdom <- inhouseDb_most_structures_2plot_wide %>%
      filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
      count(organismCleaned_dbTaxo_1kingdom) %>%
      arrange(desc(n))
    dbnumostinchi <- as.numeric(nrow(
      inhouseDbMeta %>%
        filter(structureCleanedInchikey2D == mostinchi) %>%
        distinct(database)
    ))
    upset(
      inhouseDb_most_structures_2plot_wide,
      nsets = 10,
      query.legend = "top",
      queries = list(
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_1kingdom",
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
            "organismCleaned_dbTaxo_1kingdom",
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
            "organismCleaned_dbTaxo_1kingdom",
            c(mostkingdom[1, 1])
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
      mainbar.y.label = "Unique sources per intersection",
      sets.x.label = "Unique source per database",
      set_size.show = TRUE
    )
    dev.off()
  })

  try({
    cat("drawing upset plot of quercetin repartition in sub-DBs \n")
    pdf(
      file = file.path(pathDataProcessedFigures, "quercetin.pdf"),
      width = 16,
      height = 9
    )
    mostinchi2 <- as.character(inhouseDb_most_structures[2, 1])
    inhouseDb_most_structures_2plot <-
      inhouseDbMeta %>%
      filter(structureCleanedInchikey2D == mostinchi2) %>%
      filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        database,
        .keep_all = TRUE
      ) %>%
      group_by(database, organismCleaned) %>%
      count(structureCleanedInchikey2D) %>%
      ungroup()
    inhouseDb_most_structures_2plot_wide <-
      inhouseDb_most_structures_2plot %>%
      pivot_wider(
        names_from = database,
        values_from = n
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = is.na(.),
          values = 0
        )
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = . >= 1,
          values = 1
        )
      ) %>%
      data.frame()
    inhouseDb_most_structures_2plot_wide <-
      left_join(
        inhouseDb_most_structures_2plot_wide,
        bio
      ) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        .keep_all = TRUE
      )
    mostkingdom <- inhouseDb_most_structures_2plot_wide %>%
      filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
      count(organismCleaned_dbTaxo_1kingdom) %>%
      arrange(desc(n))
    dbnumostinchi <- as.numeric(nrow(
      inhouseDbMeta %>%
        filter(structureCleanedInchikey2D == mostinchi2) %>%
        distinct(database)
    ))
    upset(
      inhouseDb_most_structures_2plot_wide,
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
      set_size.show = TRUE
    )
    dev.off()
  })

  cat("getting list of most studied plants for exploration \n")
  inhouseDb_most_plant <- inhouseDbMeta %>%
    filter(!is.na(organismCleaned)) %>%
    filter(
      organismCleaned_dbTaxo_1kingdom == "Plantae" |
        organismCleaned_dbTaxo_1kingdom == "Viridiplantae"
    ) %>%
    count(organismCleaned) %>%
    arrange(desc(n))

  cat("defining drawing function for most studied plants ... \n")
  getGraphStudiedPlant <- function(plant) {
    try({
      mostplant <-
        as.character(inhouseDb_most_plant[inhouseDb_most_plant$organismCleaned == plant, 1])
      inhouseDb_most_organism_2plot <-
        inhouseDbMeta %>%
        filter(organismCleaned == mostplant) %>%
        distinct(structureCleanedInchikey2D,
          organismCleaned,
          database,
          .keep_all = TRUE
        ) %>%
        group_by(organismCleaned, database) %>%
        count(structureCleanedInchikey2D) %>%
        ungroup()
      inhouseDb_most_organism_2plot_wide <-
        inhouseDb_most_organism_2plot %>%
        pivot_wider(
          names_from = database,
          values_from = n
        ) %>%
        mutate_at(
          .vars = c(3:ncol(.)),
          ~ replace(
            x = .,
            list = is.na(.),
            values = 0
          )
        ) %>%
        mutate_at(
          .vars = c(3:ncol(.)),
          ~ replace(
            x = .,
            list = . >= 1,
            values = 1
          )
        ) %>%
        data.frame()
      inhouseDb_most_organism_2plot_wide <-
        left_join(
          inhouseDb_most_organism_2plot_wide,
          chemo
        ) %>%
        distinct(structureCleanedInchikey2D,
          organismCleaned,
          .keep_all = TRUE
        )
      dbnumostorganism <- as.numeric(nrow(
        inhouseDbMeta %>%
          filter(organismCleaned == mostplant) %>%
          distinct(database)
      ))
      mostsuperclasses <- inhouseDb_most_organism_2plot_wide %>%
        filter(!is.na(structureCleaned_2superclass)) %>%
        count(structureCleaned_2superclass) %>%
        arrange(desc(n)) %>%
        head(10)
      upset(
        inhouseDb_most_organism_2plot_wide,
        nsets = 10,
        query.legend = "top",
        queries = list(
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1],
                mostsuperclasses[6, 1],
                mostsuperclasses[7, 1],
                mostsuperclasses[8, 1],
                mostsuperclasses[9, 1],
                mostsuperclasses[10, 1]
              )
            ),
            active = TRUE,
            color = "#6a3d9a",
            query.name = mostsuperclasses[10, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1],
                mostsuperclasses[6, 1],
                mostsuperclasses[7, 1],
                mostsuperclasses[8, 1],
                mostsuperclasses[9, 1]
              )
            ),
            active = TRUE,
            color = "#cab2d6",
            query.name = mostsuperclasses[9, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1],
                mostsuperclasses[6, 1],
                mostsuperclasses[7, 1],
                mostsuperclasses[8, 1]
              )
            ),
            active = TRUE,
            color = "#ff7f00",
            query.name = mostsuperclasses[8, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1],
                mostsuperclasses[6, 1],
                mostsuperclasses[7, 1]
              )
            ),
            active = TRUE,
            color = "#fdbf6f",
            query.name = mostsuperclasses[7, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1],
                mostsuperclasses[6, 1]
              )
            ),
            active = TRUE,
            color = "#e31a1c",
            query.name = mostsuperclasses[6, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1],
                mostsuperclasses[5, 1]
              )
            ),

            active = TRUE,
            color = "#fb9a99",
            query.name = mostsuperclasses[5, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1],
                mostsuperclasses[4, 1]
              )
            ),
            active = TRUE,
            color = "#33a02c",
            query.name = mostsuperclasses[4, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1],
                mostsuperclasses[3, 1]
              )
            ),
            active = TRUE,
            color = "#b2df8a",
            query.name = mostsuperclasses[3, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(
                mostsuperclasses[1, 1],
                mostsuperclasses[2, 1]
              )
            ),
            active = TRUE,
            color = "#1f78b4",
            query.name = mostsuperclasses[2, 1]
          ),
          list(
            query = elements,
            params = list(
              "structureCleaned_2superclass",
              c(mostsuperclasses[1, 1])
            ),
            active = TRUE,
            color = "#a6cee3",
            query.name = mostsuperclasses[1, 1]
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
        mainbar.y.label = "Unique structures per intersection",
        sets.x.label = "Unique structures per database",
        set_size.show = TRUE
      )
    })
  }

  cat("... drawing Cannabis sativa metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "cannabisSativa.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Cannabis sativa")
  dev.off()

  cat("... drawing Cryptomeria japonica metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "cryptomeriaJaponica.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Cryptomeria japonica")
  dev.off()

  cat("... drawing Tripterygium wilfordii metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "tripterygiumWilfordii.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Tripterygium wilfordii")
  dev.off()

  cat("... drawing Citrus aurantium metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "citrusAurantium.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Citrus aurantium")
  dev.off()

  cat("... drawing Camellia sinensis metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "camelliaSinensis.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Camellia sinensis")
  dev.off()

  cat("... drawing Panax ginseng metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "panaxGinseng.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Panax ginseng")
  dev.off()

  cat("... drawing Vitis vinifera metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "vitisVinifera.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Vitis vinifera")
  dev.off()

  cat("... drawing Zingiber officinale metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "zingiberOfficinale.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Zingiber officinale")
  dev.off()

  cat("... drawing Capsicum annuum metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "capsicumAnnuum.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Capsicum annuum annuum")
  dev.off()

  cat("... drawing Arabidopsis thaliana metabolites repartition \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "arabidopsisThaliana.pdf"),
    width = 16,
    height = 9
  )

  getGraphStudiedPlant(plant = "Arabidopsis thaliana")
  dev.off()

  cat("getting list of most studied chemical subclasses for exploration \n")
  inhouseDb_most_chemical_subclass <-
    inhouseDbMeta %>%
    filter(!is.na(structureCleaned_4subclass)) %>%
    count(structureCleaned_4subclass) %>%
    arrange(desc(n))

  cat("defining drawing function for most chemical subclasses ... \n")
  getGraphChemicalClass <- function(subclass) {
    inhouseDb_most_chemical_class_2plot <-
      inhouseDbMeta %>%
      filter(structureCleaned_4subclass == subclass) %>%
      distinct(structureCleanedInchikey2D,
        organismCleaned,
        database,
        .keep_all = TRUE
      ) %>%
      group_by(organismCleaned_dbTaxo_5family, database) %>%
      count(structureCleanedInchikey2D) %>%
      ungroup()

    inhouseDb_most_chemical_class_2plot_wide <-
      inhouseDb_most_chemical_class_2plot %>%
      pivot_wider(
        names_from = database,
        values_from = n
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = is.na(.),
          values = 0
        )
      ) %>%
      mutate_at(
        .vars = c(3:ncol(.)),
        ~ replace(
          x = .,
          list = . >= 1,
          values = 1
        )
      ) %>%
      data.frame()

    dbnumostchemicalclass <- as.numeric(nrow(
      inhouseDbMeta %>%
        filter(structureCleaned_4subclass == subclass) %>%
        distinct(database)
    ))

    mostfamilies <-
      inhouseDb_most_chemical_class_2plot_wide %>%
      filter(!is.na(organismCleaned_dbTaxo_5family)) %>%
      count(organismCleaned_dbTaxo_5family) %>%
      arrange(desc(n)) %>%
      head(10)

    upset(
      inhouseDb_most_chemical_class_2plot_wide,
      nsets = 10,
      query.legend = "top",
      queries = list(
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1],
              mostfamilies[6, 1],
              mostfamilies[7, 1],
              mostfamilies[8, 1],
              mostfamilies[9, 1],
              mostfamilies[10, 1]
            )
          ),
          active = TRUE,
          color = "#6a3d9a",
          query.name = mostfamilies[10, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1],
              mostfamilies[6, 1],
              mostfamilies[7, 1],
              mostfamilies[8, 1],
              mostfamilies[9, 1]
            )
          ),
          active = TRUE,
          color = "#cab2d6",
          query.name = mostfamilies[9, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1],
              mostfamilies[6, 1],
              mostfamilies[7, 1],
              mostfamilies[8, 1]
            )
          ),
          active = TRUE,
          color = "#ff7f00",
          query.name = mostfamilies[8, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1],
              mostfamilies[6, 1],
              mostfamilies[7, 1]
            )
          ),
          active = TRUE,
          color = "#fdbf6f",
          query.name = mostfamilies[7, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1],
              mostfamilies[6, 1]
            )
          ),
          active = TRUE,
          color = "#e31a1c",
          query.name = mostfamilies[6, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1],
              mostfamilies[5, 1]
            )
          ),

          active = TRUE,
          color = "#fb9a99",
          query.name = mostfamilies[5, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1],
              mostfamilies[4, 1]
            )
          ),
          active = TRUE,
          color = "#33a02c",
          query.name = mostfamilies[4, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1],
              mostfamilies[3, 1]
            )
          ),
          active = TRUE,
          color = "#b2df8a",
          query.name = mostfamilies[3, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(
              mostfamilies[1, 1],
              mostfamilies[2, 1]
            )
          ),
          active = TRUE,

          color = "#1f78b4",
          query.name = mostfamilies[2, 1]
        ),
        list(
          query = elements,
          params = list(
            "organismCleaned_dbTaxo_5family",
            c(mostfamilies[1, 1])
          ),
          active = TRUE,
          color = "#a6cee3",
          query.name = mostfamilies[1, 1]
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
      mainbar.y.label = "Unique structures per intersection",
      sets.x.label = "Unique structures per database",
      set_size.show = TRUE
    )
  }

  cat("... drawing Terpene lactones repartition \n")
  try({
    pdf(
      file = file.path(pathDataProcessedFigures, "terpeneLactones.pdf"),
      width = 16,
      height = 9
    )
    getGraphChemicalClass(subclass = "Terpene lactones")
    dev.off()
  })

  cat("... drawing Furanocoumarins repartition \n")
  try({
    pdf(
      file = file.path(pathDataProcessedFigures, "furanocoumarins.pdf"),
      width = 16,
      height = 9
    )
    getGraphChemicalClass(subclass = "Furanocoumarins")
    dev.off()
  })

  cat("drawing repartition of metabolites among biological kingdoms \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "kingdoms.pdf"),
    width = 16,
    height = 9
  )
  inhouseDb_kingdoms <- inhouseDbMeta %>%
    filter(!is.na(organismCleaned_dbTaxo_1kingdom)) %>%
    distinct(structureCleanedInchikey2D,
      organismCleaned_dbTaxo_1kingdom,
      .keep_all = TRUE
    ) %>%
    group_by(organismCleaned_dbTaxo_1kingdom) %>%
    count(structureCleanedInchikey2D) %>%
    ungroup()

  inhouseDb_kingdoms_wide <-
    inhouseDb_kingdoms %>%
    pivot_wider(
      names_from = organismCleaned_dbTaxo_1kingdom,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    data.frame()

  inhouseDb_kingdoms_wide <-
    left_join(inhouseDb_kingdoms_wide, chemo) %>%
    distinct(structureCleanedInchikey2D, .keep_all = TRUE)

  mostsuperclass <- inhouseDb_kingdoms_wide %>%
    count(structureCleaned_2superclass) %>%
    arrange(desc(n)) %>%
    head(10)

  upset(
    inhouseDb_kingdoms_wide,
    nsets = 10,
    query.legend = "top",
    queries = list(
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1],
            mostsuperclass[6, 1],
            mostsuperclass[7, 1],
            mostsuperclass[8, 1],
            mostsuperclass[9, 1],
            mostsuperclass[10, 1]
          )
        ),
        active = TRUE,
        color = "#6a3d9a",
        query.name = mostsuperclass[10, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1],
            mostsuperclass[6, 1],
            mostsuperclass[7, 1],
            mostsuperclass[8, 1],
            mostsuperclass[9, 1]
          )
        ),
        active = TRUE,
        color = "#cab2d6",
        query.name = mostsuperclass[9, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1],
            mostsuperclass[6, 1],
            mostsuperclass[7, 1],
            mostsuperclass[8, 1]
          )
        ),
        active = TRUE,
        color = "#ff7f00",
        query.name = mostsuperclass[8, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1],
            mostsuperclass[6, 1],
            mostsuperclass[7, 1]
          )
        ),
        active = TRUE,
        color = "#fdbf6f",
        query.name = mostsuperclass[7, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1],
            mostsuperclass[6, 1]
          )
        ),
        active = TRUE,
        color = "#e31a1c",
        query.name = mostsuperclass[6, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1],
            mostsuperclass[5, 1]
          )
        ),

        active = TRUE,
        color = "#fb9a99",
        query.name = mostsuperclass[5, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1],
            mostsuperclass[4, 1]
          )
        ),
        active = TRUE,
        color = "#33a02c",
        query.name = mostsuperclass[4, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1],
            mostsuperclass[3, 1]
          )
        ),
        active = TRUE,
        color = "#b2df8a",
        query.name = mostsuperclass[3, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass[1, 1],
            mostsuperclass[2, 1]
          )
        ),
        active = TRUE,
        color = "#1f78b4",
        query.name = mostsuperclass[2, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(mostsuperclass[1, 1])
        ),
        active = TRUE,
        color = "#a6cee3",
        query.name = mostsuperclass[1, 1]
      )
    ),
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 60,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per kingdom",
    set_size.show = TRUE
  )
  dev.off()

  cat("drawing repartition of metabolites among Plantae phyla \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "phyla.pdf"),
    width = 16,
    height = 9
  )

  inhouseDb_phyla <- inhouseDbMeta %>%
    filter(organismCleaned_dbTaxo_1kingdom == "Plantae") %>%
    filter(!is.na(organismCleaned_dbTaxo_2phylum)) %>%
    distinct(structureCleanedInchikey2D,
      organismCleaned_dbTaxo_2phylum,
      .keep_all = TRUE
    ) %>%
    group_by(organismCleaned_dbTaxo_2phylum) %>%
    count(structureCleanedInchikey2D) %>%
    ungroup()

  inhouseDb_phyla_wide <-
    inhouseDb_phyla %>%
    pivot_wider(
      names_from = organismCleaned_dbTaxo_2phylum,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    data.frame()

  inhouseDb_phyla_wide <-
    left_join(inhouseDb_phyla_wide, chemo) %>%
    distinct(structureCleanedInchikey2D, .keep_all = TRUE)

  mostsuperclass2 <- inhouseDb_phyla_wide %>%
    count(structureCleaned_2superclass) %>%
    arrange(desc(n)) %>%
    head(10)

  upset(
    inhouseDb_phyla_wide,
    nsets = 10,
    query.legend = "top",
    queries = list(
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1],
            mostsuperclass2[6, 1],
            mostsuperclass2[7, 1],
            mostsuperclass2[8, 1],
            mostsuperclass2[9, 1],
            mostsuperclass2[10, 1]
          )
        ),
        active = TRUE,
        color = "#6a3d9a",
        query.name = mostsuperclass2[10, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1],
            mostsuperclass2[6, 1],
            mostsuperclass2[7, 1],
            mostsuperclass2[8, 1],
            mostsuperclass2[9, 1]
          )
        ),
        active = TRUE,
        color = "#cab2d6",
        query.name = mostsuperclass2[9, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1],
            mostsuperclass2[6, 1],
            mostsuperclass2[7, 1],
            mostsuperclass2[8, 1]
          )
        ),
        active = TRUE,
        color = "#ff7f00",
        query.name = mostsuperclass2[8, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1],
            mostsuperclass2[6, 1],
            mostsuperclass2[7, 1]
          )
        ),
        active = TRUE,
        color = "#fdbf6f",
        query.name = mostsuperclass2[7, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1],
            mostsuperclass2[6, 1]
          )
        ),
        active = TRUE,
        color = "#e31a1c",
        query.name = mostsuperclass2[6, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1],
            mostsuperclass2[5, 1]
          )
        ),

        active = TRUE,
        color = "#fb9a99",
        query.name = mostsuperclass2[5, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1],
            mostsuperclass2[4, 1]
          )
        ),
        active = TRUE,
        color = "#33a02c",
        query.name = mostsuperclass2[4, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1],
            mostsuperclass2[3, 1]
          )
        ),
        active = TRUE,
        color = "#b2df8a",
        query.name = mostsuperclass2[3, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass2[1, 1],
            mostsuperclass2[2, 1]
          )
        ),
        active = TRUE,
        color = "#1f78b4",
        query.name = mostsuperclass2[2, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(mostsuperclass2[1, 1])
        ),
        active = TRUE,
        color = "#a6cee3",
        query.name = mostsuperclass2[1, 1]
      )
    ),
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 50,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per phylum",
    set_size.show = TRUE
  )
  dev.off()

  cat("drawing repartition of metabolites among Tracheophyta classes \n")
  pdf(
    file = file.path(pathDataProcessedFigures, "classes.pdf"),
    width = 16,
    height = 9
  )

  inhouseDb_classes <- inhouseDbMeta %>%
    filter(organismCleaned_dbTaxo_2phylum == "Tracheophyta") %>%
    filter(!is.na(organismCleaned_dbTaxo_3class)) %>%
    distinct(structureCleanedInchikey2D,
      organismCleaned_dbTaxo_3class,
      .keep_all = TRUE
    ) %>%
    group_by(organismCleaned_dbTaxo_3class) %>%
    count(structureCleanedInchikey2D) %>%
    ungroup()

  inhouseDb_classes_wide <-
    inhouseDb_classes %>%
    pivot_wider(
      names_from = organismCleaned_dbTaxo_3class,
      values_from = n
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = is.na(.),
        values = 0
      )
    ) %>%
    mutate_at(
      .vars = c(2:ncol(.)),
      ~ replace(
        x = .,
        list = . >= 1,
        values = 1
      )
    ) %>%
    data.frame()

  inhouseDb_classes_wide <-
    left_join(inhouseDb_classes_wide, chemo) %>%
    distinct(structureCleanedInchikey2D, .keep_all = TRUE)

  mostsuperclass3 <- inhouseDb_classes_wide %>%
    count(structureCleaned_2superclass) %>%
    arrange(desc(n)) %>%
    head(10)

  upset(
    inhouseDb_classes_wide,
    nsets = 10,
    query.legend = "top",
    queries = list(
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1],
            mostsuperclass3[6, 1],
            mostsuperclass3[7, 1],
            mostsuperclass3[8, 1],
            mostsuperclass3[9, 1],
            mostsuperclass3[10, 1]
          )
        ),
        active = TRUE,
        color = "#6a3d9a",
        query.name = mostsuperclass3[10, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1],
            mostsuperclass3[6, 1],
            mostsuperclass3[7, 1],
            mostsuperclass3[8, 1],
            mostsuperclass3[9, 1]
          )
        ),
        active = TRUE,
        color = "#cab2d6",
        query.name = mostsuperclass3[9, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1],
            mostsuperclass3[6, 1],
            mostsuperclass3[7, 1],
            mostsuperclass3[8, 1]
          )
        ),
        active = TRUE,
        color = "#ff7f00",
        query.name = mostsuperclass3[8, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1],
            mostsuperclass3[6, 1],
            mostsuperclass3[7, 1]
          )
        ),
        active = TRUE,
        color = "#fdbf6f",
        query.name = mostsuperclass3[7, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1],
            mostsuperclass3[6, 1]
          )
        ),
        active = TRUE,
        color = "#e31a1c",
        query.name = mostsuperclass3[6, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1],
            mostsuperclass3[5, 1]
          )
        ),

        active = TRUE,
        color = "#fb9a99",
        query.name = mostsuperclass3[5, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1],
            mostsuperclass3[4, 1]
          )
        ),
        active = TRUE,
        color = "#33a02c",
        query.name = mostsuperclass3[4, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1],
            mostsuperclass3[3, 1]
          )
        ),
        active = TRUE,
        color = "#b2df8a",
        query.name = mostsuperclass3[3, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(
            mostsuperclass3[1, 1],
            mostsuperclass3[2, 1]
          )
        ),
        active = TRUE,
        color = "#1f78b4",
        query.name = mostsuperclass3[2, 1]
      ),
      list(
        query = elements,
        params = list(
          "structureCleaned_2superclass",
          c(mostsuperclass3[1, 1])
        ),
        active = TRUE,
        color = "#a6cee3",
        query.name = mostsuperclass3[1, 1]
      )
    ),
    # mb.ratio = c(0.7, 0.3),
    order.by = "freq",
    nintersects = 81,
    # empty.intersections = "on",
    number.angles = 30,
    point.size = 5,
    line.size = 2,
    text.scale = 2,
    mainbar.y.label = "Unique structures per intersection",
    sets.x.label = "Unique structures per class",
    set_size.show = TRUE
  )
  dev.off()

  try({
    cat("preparing alluvial plot ... \n")
    openDbMetaValidated <- openDbMeta %>%
      mutate(validation = "validated")
    full <- left_join(openDbMaximal, openDbMetaValidated) %>%
      mutate(
        validation = ifelse(
          test = !is.na(validation),
          yes = validation,
          no = "not_validated"
        ),
        cleaned_reference = ifelse(
          test = !is.na(referenceCleanedDoi) |
            !is.na(referenceCleanedPmcid) |
            !is.na(referenceCleanedPmid),
          yes = "reference_Yes",
          no = "reference_No"
        ),
        cleaned_organism = ifelse(
          test = !is.na(organismCleaned),
          yes = "organism_Yes",
          no = "organism_No"
        ),
        cleaned_structure = ifelse(
          test = !is.na(structureCleanedSmiles) |
            !is.na(structureCleanedInchi) |
            !is.na(structureCleanedInchikey3D),
          yes = "structure_Yes",
          no = "structure_No"
        ),
      ) %>%
      distinct(
        database,
        organismOriginal,
        structureType,
        structureValue,
        referenceType,
        referenceValue,
        organismCleaned,
        structureCleanedInchikey3D,
        referenceCleanedTitle,
        cleaned_structure,
        cleaned_organism,
        cleaned_reference,
        validation,
      ) %>%
      pivot_longer(
        cols = 11:13,
        values_drop_na = TRUE
      ) %>%
      mutate(organism_originalType = "organism") %>%
      select(
        database,
        structure_originalValue = structureValue,
        structure_originalType = structureType,
        organism_originalValue = organismOriginal,
        organism_originalType,
        reference_originalValue = referenceValue,
        reference_originalType = referenceType,
        cleanedType = name,
        cleanedValue = value,
        validation
      ) %>%
      distinct() %>%
      data.frame()
    ready_1 <- full %>%
      pivot_wider(
        names_from = reference_originalType,
        names_prefix = "reference_",
        values_from = reference_originalValue,
        values_fn = first
      ) %>%
      pivot_wider(
        names_from = structure_originalType,
        names_prefix = "structure_",
        values_from = structure_originalValue,
        values_fn = first
      ) %>%
      select(
        database,
        organism_organism = organism_originalValue,
        structure_structureSmiles = structure_smiles,
        structure_structureInchi = structure_inchi,
        structure_structureNominal = structure_nominal,
        reference_referenceDoi = reference_doi,
        reference_referencePmid = reference_pubmed,
        reference_referencePublishingDetails = reference_publishingDetails,
        reference_referenceTitle = reference_title,
        reference_referenceOriginal = reference_original,
        reference_referenceSplit = reference_split,
        cleanedType,
        cleanedValue,
        validation
      )
    ready_2 <- ready_1 %>%
      pivot_longer(
        cols = 2:11,
        names_to = c("origin", "originalType"),
        names_sep = "_",
        values_to = "originalValue"
      ) %>%
      select(
        database,
        originalType,
        originalValue,
        cleanedType,
        cleanedValue,
        validation
      ) %>%
      distinct() %>%
      filter(!is.na(originalValue))
    ready_3 <- ready_2 %>%
      filter(validation == "validated") %>%
      group_by(
        database,
        validation
      ) %>%
      count(name = "count") %>%
      arrange(desc(count))
    sunk <- ready_2 %>%
      filter(database %in% ready_3$database) %>%
      group_by(
        database,
        originalType,
        cleanedType,
        validation
      ) %>%
      count(name = "count") %>%
      filter(
        gsub(
          pattern = ".*_",
          replacement = "",
          x = cleanedType
        ) %in% substr(
          x = originalType,
          start = 1,
          stop = 9
        )
      ) %>%
      ungroup() %>%
      arrange(desc(count, validation, database)) %>%
      mutate(
        validation = gsub(
          pattern = "notValidated",
          replacement = "not_validated",
          x = validation
        ),
        originalType = gsub(
          pattern = "organismOriginal",
          replacement = "organism_original",
          x = originalType
        ),
        originalType = gsub(
          pattern = "structureSmiles",
          replacement = "structure_smiles",
          x = originalType
        ),
        originalType = gsub(
          pattern = "structureInchi",
          replacement = "structure_inchi",
          x = originalType
        ),
        originalType = gsub(
          pattern = "structureNominal",
          replacement = "structure_name",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referenceDoi",
          replacement = "reference_doi",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referencePmid",
          replacement = "reference_pmid",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referencePublishingDetails",
          replacement = "reference_publishing_details",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referenceOriginal",
          replacement = "reference_original",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referenceSplit",
          replacement = "reference_split",
          x = originalType
        ),
        originalType = gsub(
          pattern = "referenceTitle",
          replacement = "reference_title",
          x = originalType
        ),
        cleanedType = gsub(
          pattern = "cleaned_",
          replacement = "",
          x = cleanedType
        ),
      )
    legend <- with(sunk, reorder(database, count))
    legend_v <- with(sunk, reorder(validation, desc(count)))
    cat("drawing alluvial \n")
    pdf(
      file = file.path(pathDataProcessedFigures, "alluvial.pdf"),
      width = 96,
      height = 54
    )
    ggplot(
      as.data.frame(sunk),
      aes(
        y = count,
        axis1 = database,
        axis2 = originalType,
        axis3 = cleanedType
      )
    ) +
      geom_stratum(decreasing = TRUE) +
      geom_alluvium(
        aes(fill = legend_v),
        aes.bind = "alluvia",
        lode.guidance = "forward",
        decreasing = TRUE
      ) +
      geom_flow(
        aes(
          fill = legend_v,
          colour = legend_v
        ),
        aes.bind = "alluvia",
        aes.flow = "forward",
        stat = after_stat("alluvium"),
        decreasing = TRUE
      ) +
      geom_fit_text(
        stat = "stratum",
        min.size = 0,
        grow = TRUE,
        width = 1 / 4,
        aes(label = after_stat(stratum)),
        decreasing = TRUE
      ) +
      scale_x_discrete(limits = c("database", "original", "cleaned")) +
      # scale_y_continuous(trans = 'log10', name = "log10(count)") +
      scale_fill_manual(values = c("#fb9a99", "#1f78b4")) +
      scale_colour_manual(values = c("#fb9a99", "#1f78b4")) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.text.y = element_text(size = rel(7)),
        axis.title.y = element_text(size = rel(1)),
        # axis.ticks = element_blank(),
        axis.text.x = element_text(size = rel(7)),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(5)),
      )
    dev.off()
  })

  detach("package:Hmisc", unload = TRUE)

  try({
    cat("drawing interactive biological tree \n")
    data4tree_bio <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_1kingdom) &
          !is.na(organismCleaned_dbTaxo_2phylum) &
          !is.na(organismCleaned_dbTaxo_3class) &
          !is.na(organismCleaned_dbTaxo_4order) &
          !is.na(organismCleaned_dbTaxo_5family) &
          !is.na(organismCleaned_dbTaxo_6genus) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      distinct(organismCleaned_dbTaxo_7species,
        .keep_all = TRUE
      ) %>%
      group_by(
        organismCleaned_dbTaxo_1kingdom,
        organismCleaned_dbTaxo_2phylum,
        organismCleaned_dbTaxo_3class,
        organismCleaned_dbTaxo_4order,
        organismCleaned_dbTaxo_5family,
        organismCleaned_dbTaxo_6genus,
        organismCleaned_dbTaxo_7species
      ) %>%
      summarize("percentage of parent" = n())
    Tree_bio <- collapsibleTreeSummary(
      df = data4tree_bio,
      hierarchy = c(
        "organismCleaned_dbTaxo_1kingdom",
        "organismCleaned_dbTaxo_2phylum",
        "organismCleaned_dbTaxo_3class",
        "organismCleaned_dbTaxo_4order",
        "organismCleaned_dbTaxo_5family",
        "organismCleaned_dbTaxo_6genus",
        "organismCleaned_dbTaxo_7species"
      ),
      root = "root",
      nodeSize = NULL,
      attribute = "percentage of parent",
      collapsed = TRUE,
      fillFun = grDevices::rainbow,
      maxPercent = 10,
      percentOfParent = TRUE,
      fontSize = 12,
      width = 4800,
      height = 9600
    )
    setwd(dir = pathDataProcessedFiguresHtml)
    htmlwidgets::saveWidget(
      widget = as_widget(Tree_bio),
      file = "Tree_bio.html"
    )
  })

  try({
    cat("drawing interactive chemical tree \n")
    data4tree_chemo <- openDbMeta %>%
      filter(!is.na(structureCleaned_5directParent)) %>%
      filter(
        !is.na(structureCleaned_1kingdom) &
          !is.na(structureCleaned_2superclass) &
          !is.na(structureCleaned_3class) &
          !is.na(structureCleaned_4subclass) &
          !is.na(structureCleaned_5directParent)
      ) %>%
      distinct(structureCleaned_5directParent, .keep_all = TRUE) %>%
      group_by(
        structureCleaned_1kingdom,
        structureCleaned_2superclass,
        structureCleaned_3class,
        structureCleaned_4subclass,
        structureCleaned_5directParent
      ) %>%
      summarize("percentage of parent" = n())
    Tree_chemo <- collapsibleTreeSummary(
      df = data4tree_chemo,
      hierarchy = c(
        "structureCleaned_1kingdom",
        "structureCleaned_2superclass",
        "structureCleaned_3class",
        "structureCleaned_4subclass",
        "structureCleaned_5directParent"
      ),
      root = "root",
      nodeSize = NULL,
      attribute = "percentage of parent",
      collapsed = TRUE,
      fillFun = grDevices::rainbow,
      maxPercent = 10,
      percentOfParent = TRUE,
      fontSize = 24,
      width = 3600,
      height = 1200
    )
    htmlwidgets::saveWidget(
      widget = as_widget(Tree_chemo),
      file = "Tree_chemo.html"
    )
  })

  cat("defining color palettes ... \n")
  cat("... small \n")
  paired_palette_sma <- c(
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928"
  )

  cat("... medium \n")
  paired_palette_med <- c(
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c"
  )

  cat("... 30 \n")
  paired_palette_30 <- c(
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c"
  )

  cat("... big \n")
  paired_palette_big <- c(
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c"
  )

  paired_palette_meg <- c(
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c"
  )

  cat("defining drawing function for chord diagrams ... \n")
  draw_chord <-
    function(data,
             biological_level,
             chemical_level,
             biological_filter_value = NULL,
             biological_filter_level = NULL,
             chemical_filter_value = NULL,
             chemical_filter_level = NULL,
             palette = paired_palette_med) {
      try({
        table <- data.frame(data)
        table <- table[!is.na(table[, biological_level]), ]
        table <- table[!is.na(table[, chemical_level]), ]
        if (!is.null(biological_filter_value)) {
          table <-
            table[table[, biological_filter_level] %in% biological_filter_value, ]
        }
        if (!is.null(chemical_filter_value)) {
          table <-
            table[table[, chemical_filter_level] %in% chemical_filter_value, ]
        }
        m1 <-
          as.data.table(table(table[, c(biological_level, chemical_level)]))
        m1 <- m1 %>%
          pivot_wider(names_from = 2, values_from = N) %>%
          unnest() %>%
          column_to_rownames(var = biological_level) %>%
          select(order(colSums(-.)))
        m2 <- t(m1) %>%
          as.data.frame() %>%
          select(order(colSums(-.)))
        m2$name <- colnames(m1)
        # colnames(m1) <- paste("chemo", colnames(m1), sep = "_")
        m1$name <- sort(colnames(m2[, 1:ncol(m2) - 1]))
        # colnames(m2)[1:ncol(m2)-1] <- paste("bio", colnames(m2[1:ncol(m2)-1]), sep = "_")
        test_3 <- full_join(m2, m1)
        test_3[is.na(test_3)] <- 0
        rownames(test_3) <- test_3$name
        test_3 <- test_3 %>% select(-name)
        test_4 <- as.matrix(test_3)
        test_5 <- test_4[colnames(test_4), colnames(test_4)]
        chord <- chorddiag(
          data = test_5,
          groupColors = palette,
          groupnamePadding = 10,
          groupThickness = 0.1,
          chordedgeColor = palette,
          groupnameFontsize = 18,
          margin = 260,
          showTooltips = FALSE,
          ticklabelFontsize = 0,
          showTicks = FALSE,
          showZeroTooltips = FALSE
        )
        return(chord)
      })
    }

  try({
    cat("... drawing big chord diagram \n")
    top_big_chord_bio <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_1kingdom) &
          !is.na(structureCleaned_2superclass)
      ) %>%
      group_by(organismCleaned_dbTaxo_1kingdom) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_1kingdom) %>%
      head(6)
    top_big_chord_bio <-
      top_big_chord_bio$organismCleaned_dbTaxo_1kingdom
    top_big_chord_chemo <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_1kingdom) &
          !is.na(structureCleaned_2superclass)
      ) %>%
      filter(organismCleaned_dbTaxo_1kingdom %in% top_big_chord_bio) %>%
      group_by(structureCleaned_2superclass) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_2superclass) %>%
      head(12)
    top_big_chord_chemo <-
      top_big_chord_chemo$structureCleaned_2superclass
    chord_big <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_1kingdom",
      chemical_level = "structureCleaned_2superclass",
      chemical_filter_value = top_big_chord_chemo,
      chemical_filter_level = "structureCleaned_2superclass",
      biological_filter_value = top_big_chord_bio,
      biological_filter_level = "organismCleaned_dbTaxo_1kingdom",
      palette = paired_palette_med
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_big),
      file = "Chord_big.html"
    )
  })

  try({
    cat("... drawing medium chord diagram \n")
    top_organism_med <- openDbMeta %>%
      filter(!is.na(organismCleaned_dbTaxo_5family)) %>%
      filter(structureCleaned_2superclass == "Alkaloids and derivatives") %>%
      group_by(organismCleaned_dbTaxo_5family) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_5family) %>%
      head(12)
    top_organism_med <-
      top_organism_med$organismCleaned_dbTaxo_5family
    top_chemo_med <- openDbMeta %>%
      filter(!is.na(organismCleaned_dbTaxo_5family)) %>%
      filter(structureCleaned_2superclass == "Alkaloids and derivatives") %>%
      filter(!is.na(structureCleaned_3class)) %>%
      filter(organismCleaned_dbTaxo_5family %in% top_organism_med) %>%
      group_by(structureCleaned_3class) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_3class) %>%
      head(18)
    top_chemo_med <- top_chemo_med$structureCleaned_3class
    chord_med <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_5family",
      chemical_level = "structureCleaned_3class",
      chemical_filter_value = top_chemo_med,
      chemical_filter_level = "structureCleaned_3class",
      biological_filter_value = top_organism_med,
      biological_filter_level = "organismCleaned_dbTaxo_5family",
      palette = paired_palette_30
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_med),
      file = "Chord_med.html"
    )
  })

  try({
    cat("... drawing small chord diagram \n")
    top_organism_sma <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_6genus == "Erythroxylum") %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(12)
    top_organism_sma <-
      top_organism_sma$organismCleaned_dbTaxo_7species
    top_chemo_sma <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_6genus == "Erythroxylum") %>%
      filter(organismCleaned_dbTaxo_7species %in% top_organism_sma) %>%
      filter(!is.na(structureCleaned_5directParent)) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent) %>%
      head(12)
    top_chemo_sma <- top_chemo_sma$structureCleaned_5directParent
    chord_sma <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_sma,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_sma,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_big
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_sma),
      file = "Chord_sma.html"
    )
  })

  try({
    cat("... drawing Ranunculaceae chord diagram \n")
    top_organism_ranunculaceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Ranunculaceae") %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(12)
    top_organism_ranunculaceae <-
      top_organism_ranunculaceae$organismCleaned_dbTaxo_7species
    top_chemo_ranunculaceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Ranunculaceae") %>%
      filter(organismCleaned_dbTaxo_7species %in% top_organism_ranunculaceae) %>%
      filter(!is.na(structureCleaned_5directParent)) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent) %>%
      head(12)
    top_chemo_ranunculaceae <-
      top_chemo_ranunculaceae$structureCleaned_5directParent
    chord_ranunculaceae <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_ranunculaceae,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_ranunculaceae,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_big
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_ranunculaceae),
      file = "Chord_ranunculaceae.html"
    )
  })

  try({
    cat("... drawing Papaveraceae chord diagram \n")
    top_organism_papaveraceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Papaveraceae") %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(12)
    top_organism_papaveraceae <-
      top_organism_papaveraceae$organismCleaned_dbTaxo_7species
    top_chemo_papaveraceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Papaveraceae") %>%
      filter(organismCleaned_dbTaxo_7species %in% top_organism_papaveraceae) %>%
      filter(!is.na(structureCleaned_5directParent)) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent) %>%
      head(12)
    top_chemo_papaveraceae <-
      top_chemo_papaveraceae$structureCleaned_5directParent
    chord_papaveraceae <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_papaveraceae,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_papaveraceae,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_big
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_papaveraceae),
      file = "Chord_papaveraceae.html"
    )
  })

  try({
    cat("... drawing Gentianaceae chord diagram \n")
    top_organism_gentianaceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Gentianaceae") %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(36)
    top_organism_gentianaceae <-
      top_organism_gentianaceae$organismCleaned_dbTaxo_7species
    top_chemo_gentianaceae <- openDbMeta %>%
      filter(
        !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      filter(organismCleaned_dbTaxo_5family == "Gentianaceae") %>%
      filter(structureCleaned_5directParent == "Xanthones") %>%
      filter(organismCleaned_dbTaxo_7species %in% top_organism_gentianaceae) %>%
      filter(!is.na(structureCleanedInchikey2D)) %>%
      group_by(structureCleanedInchikey2D) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleanedInchikey2D) %>%
      head(144)
    top_chemo_gentianaceae <-
      top_chemo_gentianaceae$structureCleanedInchikey2D
    chord_gentianaceae <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleanedInchikey2D",
      chemical_filter_value = top_chemo_gentianaceae,
      chemical_filter_level = "structureCleanedInchikey2D",
      biological_filter_value = top_organism_gentianaceae,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_big
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_gentianaceae),
      file = "Chord_gentianaceae.html"
    )
  })

  cat("... drawing top N chord diagrams ... \n")
  try({
    cat("... top 06 \n")
    top_organism_06 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(6)
    top_organism_06 <-
      top_organism_06$organismCleaned_dbTaxo_7species
    top_chemo_06 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(structureCleaned_5directParent)
      ) %>%
      filter(organismCleaned %in% top_organism_06) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent,
        .keep_all = TRUE
      ) %>%
      head(6)
    top_chemo_06 <- top_chemo_06$structureCleaned_5directParent
    chord_06 <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_06,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_06,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_sma
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_06),
      file = "Chord_06.html"
    )
  })

  try({
    cat("... top 12 \n")
    top_organism_12 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(12)
    top_organism_12 <-
      top_organism_12$organismCleaned_dbTaxo_7species
    top_chemo_12 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(structureCleaned_5directParent)
      ) %>%
      filter(organismCleaned %in% top_organism_12) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent,
        .keep_all = TRUE
      ) %>%
      head(12)
    top_chemo_12 <- top_chemo_12$structureCleaned_5directParent
    chord_12 <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_12,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_12,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_big
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_12),
      file = "Chord_12.html"
    )
  })

  try({
    cat("... top 24 \n")
    top_organism_24 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(organismCleaned_dbTaxo_7species)
      ) %>%
      group_by(organismCleaned_dbTaxo_7species) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(organismCleaned_dbTaxo_7species) %>%
      head(24)
    top_organism_24 <-
      top_organism_24$organismCleaned_dbTaxo_7species
    top_chemo_24 <- openDbMeta %>%
      filter(
        !is.na(organismCleaned_dbTaxo_7species) &
          !is.na(structureCleaned_5directParent) &
          !is.na(structureCleaned_5directParent)
      ) %>%
      filter(organismCleaned %in% top_organism_24) %>%
      group_by(structureCleaned_5directParent) %>%
      add_count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      distinct(structureCleaned_5directParent,
        .keep_all = TRUE
      ) %>%
      head(24)
    top_chemo_24 <- top_chemo_24$structureCleaned_5directParent
    chord_24 <- draw_chord(
      data = openDbMeta,
      biological_level = "organismCleaned_dbTaxo_7species",
      chemical_level = "structureCleaned_5directParent",
      chemical_filter_value = top_chemo_24,
      chemical_filter_level = "structureCleaned_5directParent",
      biological_filter_value = top_organism_24,
      biological_filter_level = "organismCleaned_dbTaxo_7species",
      palette = paired_palette_meg
    )
    htmlwidgets::saveWidget(
      widget = as_widget(chord_24),
      file = "Chord_24.html"
    )
  })
}

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")