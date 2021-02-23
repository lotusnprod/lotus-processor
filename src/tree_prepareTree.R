cat("This script prepares data for wonderful tree plotting \n")
cat("It currently needs 'tree_classyfireTaxonomy.R' to be run before. \n")

start <- Sys.time()

cat("sourcing ... \n")
cat("... paths \n")
source("paths.R")

cat("... libraries \n")
library(treeio)
library(ggtree)
library(ggplot2)
library(data.table)
library(DBI)
library(rotl)
library(RSQLite)
library(tidyverse)
library(vroom)

db_path <- pathDataInterimDictionariesOrganismDictionaryOTL

drv <- SQLite()

db <- dbConnect(
  drv = drv,
  dbname = db_path
)

names <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_names"
) %>%
  mutate_all(as.character)

otl <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_otl"
) %>%
  mutate_all(as.character)

meta <- dbGetQuery(
  conn = db,
  statement = "SELECT * FROM taxa_meta"
) %>%
  mutate_all(as.character)

pairs <-
  vroom_read_safe(path = pathDataInterimTablesAnalysedPlatinum)

metadata_bio <- left_join(names, otl) %>%
  left_join(., meta, by = c("ott_id" = "id")) %>%
  filter(
    rank %in% c(
      "domain",
      "kingdom",
      "phylum",
      "class",
      "order",
      "infraorder",
      "family",
      "subfamily",
      "tribe",
      "subtribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  ) %>%
  distinct() %>%
  map_df(rev) %>% ## feeling it is better that way
  distinct(ott_id, rank, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = "rank",
    values_from = c("name", "unique_name.y", "ott_id.y")
  ) %>%
  select(
    canonical_name,
    ott_id,
    Domain = name_domain,
    Kingdom = name_kingdom,
    Phylum = name_phylum,
    Class = name_class,
    Order = name_order,
    Infraorder = name_infraorder,
    Family = name_family,
    Subfamily = name_subfamily,
    Tribe = name_tribe,
    Subtribe = name_subtribe,
    Genus = name_genus,
    Subgenus = name_subgenus,
    Species = name_species,
    Subspecies = name_subspecies
  ) %>%
  map_df(rev) %>%
  coalesce()

metadata_che <- classy_temp

metadata_che_new <-
  fread(
    file = "~/GitLab/lotus/lotusProcessor/data/interim/dictionaries/structure/npclassifier/smiles_np_classified.tsv.gz",
    fill = TRUE,
    na.strings = c("", NA),
    strip.white = TRUE,
    colClasses = cols(.default = "c")
  ) %>%
  filter(!is.na(pathway)) %>%
  filter(
    pathway != "Fatty acids" |
      superclass != "Fatty acyls" |
      class != "Halogenated hydrocarbons" |
      grepl(pattern = "Cl|Br|I|F", x = smiles) ## because actually returned instead of NA/null
  ) %>%
  distinct()

pairs_metadata <- pairs %>%
  distinct(
    organismCleaned,
    structureCleanedSmiles,
    structureCleanedInchikey,
    structureCleaned_inchikey2D
  ) %>%
  left_join(., metadata_bio, by = c("organismCleaned" = "canonical_name")) %>%
  left_join(., metadata_che, by = c("structureCleanedInchikey" = "inchikey")) %>%
  left_join(., metadata_che_new, by = c("structureCleanedSmiles" = "smiles")) %>%
  distinct(organismCleaned,
    structureCleaned_inchikey2D,
    .keep_all = TRUE
  )

families <- pairs_metadata %>%
  filter(!is.na(Family)) %>%
  distinct(Family)

families_matched <-
  tnrs_match_names(names = families$Family, do_approximate_matching = FALSE)

# ott_in_tree <-
#   ott_id(families_matched)[is_in_tree(ott_id(families_matched))]

# tr <- tol_induced_subtree(ott_ids = ott_in_tree)

sitosterol <- pairs_metadata %>%
  mutate(
    structureCleaned_inchikey2D = ifelse(
      test = structureCleaned_inchikey2D == "KZJWDPNRJALLNS",
      yes = structureCleaned_inchikey2D,
      no = NA
    )
  ) %>%
  group_by(Family) %>%
  fill(structureCleaned_inchikey2D, .direction = "downup") %>%
  distinct(Family, structureCleaned_inchikey2D)

stigmastanes <- pairs_metadata %>%
  mutate(
    `04subclass` = ifelse(
      test = `04subclass` == "Stigmastanes and derivatives",
      yes = `04subclass`,
      no = NA
    )
  ) %>%
  group_by(Family) %>%
  fill(`04subclass`, .direction = "downup") %>%
  distinct(
    Family,
    `04subclass`
  )

steroids <- pairs_metadata %>%
  mutate(
    `03class` = ifelse(
      test = `03class` == "Steroids and steroid derivatives",
      yes = `03class`,
      no = NA
    )
  ) %>%
  group_by(Family) %>%
  fill(`03class`, .direction = "downup") %>%
  distinct(
    Family,
    `03class`
  )

lipids <- pairs_metadata %>%
  mutate(
    `02superclass` = ifelse(
      test = `02superclass` == "Lipids and lipid-like molecules",
      yes = `02superclass`,
      no = NA
    )
  ) %>%
  group_by(Family) %>%
  fill(`02superclass`, .direction = "downup") %>%
  distinct(
    Family,
    `02superclass`
  )

families_matched <- families_matched %>%
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

sitosterol_matched <-
  left_join(families_matched,
    sitosterol,
    by = c("unique_name" = "Family")
  ) %>%
  mutate(
    structureCleaned_inchikey2D = ifelse(
      test = !is.na(structureCleaned_inchikey2D),
      yes = structureCleaned_inchikey2D,
      no = "noSitosterol"
    )
  ) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

families_restricted <- pairs_metadata %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  add_count() %>%
  ungroup() %>%
  filter(n >= 100) %>%
  distinct(Family)

families_matched_restricted <-
  tnrs_match_names(
    names = families_restricted$Family,
    do_approximate_matching = FALSE
  )

ott_in_tree_restricted <-
  ott_id(families_matched_restricted)[is_in_tree(ott_id(families_matched_restricted))]

families_restricted <- families_restricted %>%
  filter(Family %in% names(ott_in_tree_restricted))

families_matched_restricted <-
  tnrs_match_names(names = families_restricted$Family)

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

tr_restricted <-
  tol_induced_subtree(ott_ids = ott_in_tree_restricted)

sitosterol_matched_restricted <-
  left_join(families_matched_restricted,
    sitosterol,
    by = c("unique_name" = "Family")
  ) %>%
  mutate(
    structureCleaned_inchikey2D = ifelse(
      test = !is.na(structureCleaned_inchikey2D),
      yes = structureCleaned_inchikey2D,
      no = "noSitosterol"
    )
  ) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

stigmastanes_matched <-
  left_join(families_matched_restricted,
    stigmastanes,
    by = c("unique_name" = "Family")
  ) %>%
  mutate(`04subclass` = ifelse(
    test = !is.na(`04subclass`),
    yes = `04subclass`,
    no = "znoStigmastanes"
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

steroids_matched <-
  left_join(families_matched_restricted,
    steroids,
    by = c("unique_name" = "Family")
  ) %>%
  mutate(`03class` = ifelse(
    test = !is.na(`03class`),
    yes = `03class`,
    no = "znoSteroids"
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

lipids_matched <-
  left_join(families_matched_restricted,
    lipids,
    by = c("unique_name" = "Family")
  ) %>%
  mutate(`02superclass` = ifelse(
    test = !is.na(`02superclass`),
    yes = `02superclass`,
    no = "znoLipids"
  )) %>%
  relocate(key, .before = search_string) %>%
  data.frame()

# ggtree(tr, layout = "circular")
#
# ggtree(tr, layout = "circular") %<+% sitosterol_matched +
#   geom_tree(aes(color = structureCleaned_inchikey2D))

ggtree(tr_restricted, layout = "circular") %<+% sitosterol_matched_restricted +
  geom_tree(aes(color = structureCleaned_inchikey2D))

ggtree(tr = tr_restricted, layout = "circular") %<+% stigmastanes_matched +
  aes(color = X04subclass)

ggtree(tr = tr_restricted, layout = "circular") %<+% steroids_matched +
  aes(color = X03class)

ggtree(tr = tr_restricted, layout = "circular") %<+% lipids_matched +
  aes(color = X02superclass)

end <- Sys.time()

cat("Script finished in", format(end - start), "\n")
