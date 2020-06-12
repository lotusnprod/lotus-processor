#title: "Chemo hierarchizeR"

#loading functions
source("functions.R")

#writing paths
##inputs
###iteration name
iteration_name <- "first_run"

###iteration num
iteration_num <- "4"

###db chemo
inpath_chemo <- "outputs/tables/3_curated/curatedStructure.tsv.zip"

###classyfire input
inpath_classy_in <- paste("Z_old_to_keep/XX_OLD_06_classyfire/input/",
                 iteration_name,
                 "_",
                 iteration_num,
                 ".tsv",
                 sep = "")

###classyfire output(s)
inpath_classy <- paste("Z_old_to_keep/XX_OLD_06_classyfire/output/",
                       iteration_name,
                       "_classy_",
                       iteration_num,
                       ".tsv",
                       sep = "")

###already classified InChIs
inpath_inchi <- "outputs/tables/inchi.tsv"

##output
outpath <- paste("outputs/tables/4_additional_metadata/classyfire/",
                 iteration_name,
                 "_classy_",
                 iteration_num,
                 ".tsv",
                 sep = "")

#loading files
##classyfire input
classy_in <- read_delim(
  file = inpath_classy_in,
  delim = "\t",
  col_names = FALSE,
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  mutate_all(as.character) %>%
  select(CompoundID = 1,
         inchi = 2)

##db chemo sanitized
data_chemo_sanitized <- read_delim(
  file = gzfile(inpath_chemo),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

##classyfire output(s)
data_classy <- read_delim(
  file = inpath_classy,
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)%>%
  mutate_all(as.character) %>%
  select(-ParentName) %>%
  group_by(CompoundID) %>%
  mutate(m = seq_along(CompoundID))

##already classified InChIs
inchi <- read_delim(
  file = inpath_inchi,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) %>% 
  mutate_all(as.character)

##classyfire ontologies
ontologies <- as.data.frame(fromJSON(txt = "Z_old_to_keep/XX_OLD_06_classyfire/tax_nodes.json", simplifyDataFrame = TRUE))

#regenerating proper classyfire ontologies
ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id))

ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id))

ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id))

ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id))

ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id))

ontologies <-
  full_join(ontologies,
            ontologies,
            by = c("chemont_id" = "parent_chemont_id")) %>%
  select(-chemont_id, chemont_id = chemont_id.y) %>%
  filter(!is.na(chemont_id)) %>%
  select(
    ChemOntID = 66,
    "01" = 54,
    "02" = 55,
    "03" = 56,
    "04" = 57,
    "05" = 58,
    "06" = 59,
    "07" = 60,
    "08" = 61,
    "09" = 62,
    "10" = 63,
    "11" = 64,
    "12" = 65
  )

ontologies$new <-
  apply((ontologies), 1, function(x)
    paste(x[!is.na(x)], collapse = "|"))

ontologies <- ontologies %>%
  select(Classyfire = new) %>%
  cSplit(splitCols = "Classyfire",
         sep = "|")

colnames(ontologies)[1] <- "ChemOntID"

ontologies$ParentName <-
  apply(ontologies, 1, function(x)
    tail(na.omit(x), 1))

ontologies <- ontologies %>%
  select(
    ChemOntID,
    chemo_lowertaxon = 14,
    chemo_01_kingdom = 3,
    chemo_02_superclass = 4,
    chemo_03_class = 5,
    chemo_04_subclass = 6,
    chemo_05_parent = 7,
    chemo_06_subparent_1 = 8,
    chemo_07_subparent_2 = 9,
    chemo_08_subparent_3 = 10,
    chemo_09_subparent_4 = 11,
    chemo_10_subparent_5 = 12,
    chemo_11_subparent_6 = 13,
  )

#joining to classyfire output
data_classy_onto <- left_join(data_classy, ontologies) %>% 
  select(chemo_chemontid = ChemOntID,
         everything())

#joining to data chemo
data_classy_final <- left_join(classy_in, data_classy_onto) %>%
  filter(m == 1) %>%
  select(-m, -CompoundID) %>%
  filter(!is.na(inchi)) %>%
  mutate_all(as.character)

#adding already classyfied InChI's
data_chemo_final <- full_join(data_classy_final, inchi)

#adding rdkit terms
data_chemo_final <- left_join(data_chemo_sanitized,
                              data_chemo_final) 

#exporting
write.table(
  x = data_chemo_final,
  file = outpath,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
