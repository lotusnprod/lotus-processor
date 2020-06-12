#title: "PAMDB scrapeR"

##functions
source("../../functions.R")

outpath <- "0_initial_files/PAMDB_scraped.tsv.zip"

url <- 'http://pseudomonas.umaryland.edu/PAMDB?MetID=PAMDB'

X <- (c(1:7100, 100000:100999, 110000:110999, 120000:120999))

getpamdb <- function(X)
{
  tryCatch({
    cd_id <- str_pad (X, 6, pad = "0")
    url_id <- paste(url, cd_id)
    url_id <- gsub("\\s", "", url_id)
    df1 <- read_html(url_id) %>%
      html_node(xpath = "body/main/table") %>%
      html_table(., fill = TRUE)
    df2 <- tibble(df1[1:17, 1], df1[1:17, 2])
  },
  error = function(e) {
    "Timed out!"
  })
}

PAMDB <- invisible(
  pbmclapply(
    FUN = getpamdb,
    X = X,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE,
    mc.silent = TRUE,
    mc.cores = (parallel::detectCores() - 2),
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
)

PAMDB_2 <- PAMDB[PAMDB != "Timed out!"]

PAMDB_3 <- bind_rows(PAMDB_2)

PAMDB_3$level <- as.numeric(gl(nrow(PAMDB_3) / 17, 17))

colnames(PAMDB_3) <- c("name", "value", "level")

PAMDB_3$name <- y_as_na(PAMDB_3$name, "")
PAMDB_3$value <- y_as_na(PAMDB_3$value, "")

PAMDB_4 <- PAMDB_3 %>%
  filter(!is.na(name)) %>%
  group_by(level) %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  ungroup() %>%
  select(
    Version,
    UpdateDate = `Update Date`,
    MetaboliteID = `Metabolite ID`,
    Identification,
    Name = `Name:`,
    Description = `Description:`,
    Structure,
    Synonyms = `Synonyms:`,
    ChemicalFormula = `Chemical Formula:`,
    AverageMolecularWeight = `Average Molecular Weight:`,
    MonoisotopicMolecularWeight = 12,
    InChIKey = `InChI Key:`,
    InChI = `InChI:`,
    CAS_number = 15,
    IUPAC_Name = `IUPAC Name:`,
    Traditional_Name = `Traditional IUPAC Name:`,
    SMILES = `SMILES:`,
    ChemicalFormula_2 = `Chemical Formula`,
    AverageMolecularWeight_2 = `Average Molecular Weight`,
    MonoisotopicMolecularWeight_2 = `Monoisotopic Molecular Weight`,
    IUPAC_Name_2 = `IUPAC Name`,
    Traditional_Name_2 = `Traditional Name`,
    CAS_number_2 = `CAS Registry Number`,
    SMILES_2 = SMILES,
    InChI_2 = `InChI Identifier`,
    InChIKey_2 = `InChI Key`
  )

PAMDB_5 <- PAMDB_4 %>%
  rowwise() %>%
  mutate(
    ChemicalFormula = ifelse(is.na(ChemicalFormula), ChemicalFormula_2, ChemicalFormula),
    AverageMolecularWeight = ifelse(
      is.na(AverageMolecularWeight),
      AverageMolecularWeight_2,
      AverageMolecularWeight
    ),
    MonoisotopicMolecularWeight = ifelse(
      is.na(MonoisotopicMolecularWeight),
      MonoisotopicMolecularWeight_2,
      MonoisotopicMolecularWeight
    ),
    IUPAC_Name = ifelse(is.na(IUPAC_Name), IUPAC_Name_2, IUPAC_Name),
    Traditional_Name = ifelse(
      is.na(Traditional_Name),
      Traditional_Name_2,
      Traditional_Name
    ),
    CAS_number = ifelse(is.na(CAS_number), CAS_number_2, CAS_number),
    SMILES = ifelse(is.na(SMILES), SMILES_2, SMILES),
    InChI = ifelse(is.na(InChI), InChI_2, InChI),
    InChIKey = ifelse(is.na(InChIKey), InChIKey_2, InChIKey)
  ) %>%
  ungroup() %>%
  select(
    Version,
    UpdateDate,
    MetaboliteID,
    Identification,
    Name,
    Description,
    Structure,
    Synonyms,
    ChemicalFormula,
    AverageMolecularWeight,
    MonoisotopicMolecularWeight,
    InChIKey,
    InChI,
    CAS_number,
    IUPAC_Name,
    Traditional_Name,
    SMILES
  )

PAMDB_5[] <- lapply(PAMDB_5, function(x)
  gsub("\r\n", " ", x))
PAMDB_5[] <- lapply(PAMDB_5, function(x)
  gsub("\r", " ", x))
PAMDB_5[] <- lapply(PAMDB_5, function(x)
  gsub("\n", " ", x))

#exporting
write.table(
  x = PAMDB_5,
  file = gzfile(description = outpath,
                compression = 9,
                encoding = "UTF-8"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  fileEncoding = "UTF-8"
)
