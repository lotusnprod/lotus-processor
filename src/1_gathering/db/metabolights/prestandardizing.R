# title: "Metabolights precleaneR"

# loading paths
source("paths.R")

# loading functions
source("functions.R")

data <- xmlParse(pathDataExternalDbSourceMetabolightsComplete)

xml_data <- xmlToList(data)

x <- 1:length(xml_data[["entries"]])

getid <- function(x) {
  i <- xml_data[["entries"]][[x]][[".attrs"]][["id"]]
  return(i)
}

id <- pbmclapply(
  FUN = getid,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

getadditionalfields <- function(x) {
  j <- xml_data[["entries"]][[x]][["additional_fields"]]
  return(j)
}

additionalfields <- pbmclapply(
  FUN = getadditionalfields,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

getref <- function(x) {
  h <- xml_data[["entries"]][[x]][["cross_references"]][["ref"]]
  return(h)
}

references <- pbmclapply(
  FUN = getref,
  X = x,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = TRUE,
  mc.cores = (parallel::detectCores() - 2),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE
)

data <- tibble(id, additionalfields, references) %>%
  filter(grepl("MTBLC", id)) %>%
  unnest_wider(references) %>%
  select(
    uniqueid = id,
    refkey = dbkey,
    refdb = dbname,
    additionalfields
  )

data_clean <- data  %>%
  unnest(additionalfields) %>%
  unnest(additionalfields) %>%
  mutate(value = lag(additionalfields, 1))

data_clean$uniqueid <-
  as.character(data_clean$uniqueid)

data_clean$additionalfields <-
  as.character(data_clean$additionalfields)

data_clean_2 <- data_clean[-1, ] %>%
  filter(row_number() %% 2 == 1) %>%
  group_by(uniqueid, additionalfields) %>%
  add_count() %>%
  ungroup() %>%
  mutate(name = paste(additionalfields, n, sep = "_")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  data.frame()

data_clean_inchi <- data_clean_2 %>%
  filter(additionalfields == "inchi") %>%
  select(uniqueid, refkey, refdb, inchi_1) %>%
  unnest(4)

data_clean_iupac <- data_clean_2 %>%
  filter(additionalfields == "iupac") %>%
  select(uniqueid, iupac_1) %>%
  unnest(2)

data_clean_organism_01 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_1) %>%
  unnest(2)

data_clean_organism_02 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_2) %>%
  unnest(2)

data_clean_organism_03 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_3) %>%
  unnest(2)

data_clean_organism_04 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_4) %>%
  unnest(2)

data_clean_organism_05 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_5) %>%
  unnest(2)

data_clean_organism_06 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_6) %>%
  unnest(2)

data_clean_organism_07 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_7) %>%
  unnest(2)

data_clean_organism_08 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_8) %>%
  unnest(2)

data_clean_organism_09 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_9) %>%
  unnest(2)

data_clean_organism_10 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_10) %>%
  unnest(2)

data_clean_organism_11 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_11) %>%
  unnest(2)

data_clean_organism_12 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_12) %>%
  unnest(2)

# data_clean_organism_13 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_13) %>%
#   unnest(2)

data_clean_organism_14 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_14) %>%
  unnest(2)

data_clean_organism_15 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_15) %>%
  unnest(2)

data_clean_organism_16 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_16) %>%
  unnest(2)

data_clean_organism_17 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_17) %>%
  unnest(2)

data_clean_organism_18 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_18) %>%
  unnest(2)

# data_clean_organism_19 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_19) %>%
#   unnest(2)

data_clean_organism_20 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_20) %>%
  unnest(2)

# data_clean_organism_21 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_21) %>%
#   unnest(2)

data_clean_organism_22 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_22) %>%
  unnest(2)

# data_clean_organism_23 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_23) %>%
#   unnest(2)

data_clean_organism_24 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_24) %>%
  unnest(2)

# data_clean_organism_25 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_25) %>%
#   unnest(2)
#
# data_clean_organism_26 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_26) %>%
#   unnest(2)

# data_clean_organism_27 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_27) %>%
#   unnest(2)

# data_clean_organism_28 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_28) %>%
#   unnest(2)

# data_clean_organism_29 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_29) %>%
#   unnest(2)

# data_clean_organism_30 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_30) %>%
#   unnest(2)

# data_clean_organism_31 <- data_clean_2 %>%
#   filter(additionalfields == "organism") %>%
#   select(uniqueid, organism_31) %>%
#   unnest(2)

data_clean_organism_32 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_32) %>%
  unnest(2)

data_clean_organism_47 <- data_clean_2 %>%
  filter(additionalfields == "organism") %>%
  select(uniqueid, organism_47) %>%
  unnest(2)

data_joined <- left_join(data_clean_inchi, data_clean_iupac)
data_joined <- left_join(data_joined, data_clean_organism_01)
data_joined <- left_join(data_joined, data_clean_organism_02)
data_joined <- left_join(data_joined, data_clean_organism_03)
data_joined <- left_join(data_joined, data_clean_organism_04)
data_joined <- left_join(data_joined, data_clean_organism_05)
data_joined <- left_join(data_joined, data_clean_organism_06)
data_joined <- left_join(data_joined, data_clean_organism_07)
data_joined <- left_join(data_joined, data_clean_organism_08)
data_joined <- left_join(data_joined, data_clean_organism_09)
data_joined <- left_join(data_joined, data_clean_organism_10)
data_joined <- left_join(data_joined, data_clean_organism_11)
data_joined <- left_join(data_joined, data_clean_organism_12)
data_joined <- left_join(data_joined, data_clean_organism_14)
data_joined <- left_join(data_joined, data_clean_organism_15)
data_joined <- left_join(data_joined, data_clean_organism_16)
data_joined <- left_join(data_joined, data_clean_organism_17)
data_joined <- left_join(data_joined, data_clean_organism_18)
data_joined <- left_join(data_joined, data_clean_organism_20)
data_joined <- left_join(data_joined, data_clean_organism_22)
data_joined <- left_join(data_joined, data_clean_organism_24)
data_joined <- left_join(data_joined, data_clean_organism_32)
data_joined <- left_join(data_joined, data_clean_organism_47)

data_clean_3 <- data_joined %>%
  pivot_longer(6:ncol(.)) %>%
  filter(value != "NULL")

data_clean_final <- data_clean_3 %>%
  select(
    uniqueid,
    inchi = inchi_1,
    name = iupac_1,
    biologicalsource = value,
    reference = refkey
  )

data_clean_final[] <-
  lapply(data_clean_final, function(x)
    gsub("\r\n", " ", x))
data_clean_final[] <-
  lapply(data_clean_final, function(x)
    gsub("\r", " ", x))
data_clean_final[] <-
  lapply(data_clean_final, function(x)
    gsub("\n", " ", x))
data_clean_final[] <-
  lapply(data_clean_final, function(x)
    gsub("\t", " ", x))

write.table(
  x = data_clean_final,
  file = gzfile(
    description = pathDataExternalDbSourceMetabolightsPrecleaned,
    compression = 9,
    encoding = "UTF-8"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
