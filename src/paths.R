###############################################################################
##################################   SSOT   ###################################
###############################################################################

ssot_access <- FALSE

dbname <- "lotus"
user <- "lotusadmin"
host <- "10.9.0.1"
port <- 5432
password <- Sys.getenv("LOTUS_DB_PWD")

###############################################################################
##################################   Modes   ##################################
###############################################################################

mode <- Sys.getenv("MODE", unset = "test")
if (exists("mode_test")) {
  mode <- "test"
}
if (exists("mode_custom")) {
  mode <- "custom"
}

works_locally_only <- TRUE
safety <- FALSE

locales <- readr::locale(
  date_names = "en",
  date_format = "%AD",
  time_format = "%AT",
  decimal_mark = ".",
  grouping_mark = ",",
  tz = "UTC",
  encoding = "UTF-8",
  asciify = FALSE
)

###############################################################################
##################################   Paths   ##################################
###############################################################################
source("r/database.R")
source("r/create_dir.R")

## for steps where paths need to be adapted locally or related to other repos
molconvertPath <- "~/../../Applications/MarvinSuite/bin/molconvert"

## can be easily replaced locally
wikidataLotusExporterDataOutputPath <-
  "../../lotus-wikidata-interact/downloadLotus/data/output"
# "https://zenodo.org/record/5793224/files"

wikidataLotusExporterDataOutputTaxaPath <-
  file.path(wikidataLotusExporterDataOutputPath, "taxa.tsv")

wikidataLotusExporterDataOutputStructuresPath <-
  file.path(wikidataLotusExporterDataOutputPath, "compounds.tsv")

wikidataLotusExporterDataOutputReferencesPath <-
  file.path(wikidataLotusExporterDataOutputPath, "references.tsv")

wikidataLotusExporterDataOutputTriplesPath <-
  file.path(wikidataLotusExporterDataOutputPath, "compound_reference_taxon.tsv")

# databases for which we have no right to disseminate the data
forbidden_export <- c("antibase", "antimarin", "dnp", "foodb")

# root
## bin
pathBin <- Sys.getenv("BIN_PATH",
  unset = "../bin"
)

## opsin
pathBinOpsin <- file.path(
  pathBin,
  "opsin-cli-2.8.0-jar-with-dependencies.jar"
)

## data
pathData <- Sys.getenv("DATA_PATH",
  unset = "../data"
)

### external
pathDataExternal <-
  file.path(
    pathData,
    "external"
  )

#### db source
pathDataExternalDbSource <-
  file.path(
    pathDataExternal,
    "dbSource"
  )

### interim
pathDataInterim <-
  file.path(
    pathData,
    "interim"
  )

### processed
pathDataProcessed <-
  file.path(
    pathData,
    "processed"
  )

#### db
if (mode != "custom") {
  pathDataInterimDb <- file.path(
    pathDataInterim,
    "db"
  )
} else {
  pathDataInterimDb <- file.path(pathDataInterim, "custom")
}

databases <-
  Databases$new(
    pathDbSource = pathDataExternalDbSource,
    pathDbInterim = pathDataInterimDb
  )

databases$add(
  name = "afrotryp",
  sourceFiles = list(tsv = "afrotryp.tsv.zip"),
  interimFile = "afrotryp.tsv.gz"
)

databases$add(
  name = "alkamid",
  sourceFiles = list(
    tsv = "alkamidScraped.tsv.gz",
    tsvRef = "alkamidRefScraped.tsv.gz"
  ),
  interimFile = "alkamid.tsv.gz"
)

databases$add(
  name = "anpdb",
  sourceFiles = list(tsv = "anpdbScraped.tsv.gz"),
  interimFile = "anpdb.tsv.gz"
)

databases$add(
  name = "antibase",
  sourceFiles = list(
    sdf = "ANTIBASE_2012_FORM2.sdf",
    smi = "antibaseConverted.smi",
    tsv = "antibaseConverted.tsv.gz"
  ),
  interimFile = "antibase.tsv.gz"
)

databases$add(
  name = "antimarin",
  sourceFiles = list(
    sdf = "antimarin0311_test.sdf",
    smi = "antimarinConverted.smi",
    tsv = "antimarinConverted.tsv.gz"
  ),
  interimFile = "antimarin.tsv.gz"
)

databases$add(
  name = "biofacquim",
  sourceFiles = list(
    sdf = "BIOFACQUIM_V2.sdf",
    tsv = "biofacquim.tsv.gz"
  ),
  interimFile = "biofacquim.tsv.gz"
)

databases$add(
  name = "biophytmol",
  sourceFiles = list(tsv = "biophytmolScraped.tsv.gz"),
  interimFile = "biophytmol.tsv.gz"
)

databases$add(
  name = "carotenoiddb",
  sourceFiles = list(
    tsv = "carotenoiddbScraped.tsv.gz",
    tsvInchi = "Carotenoids_InChI_InChIKey.tsv"
  ),
  interimFile = "carotenoiddb.tsv.gz"
)

databases$add(
  name = "coconut",
  sourceFiles = list(
    sdf = "COCONUT_DB.sdf",
    tsv = "coconutConverted.tsv.gz"
  ),
  interimFile = "coconut.tsv.gz"
)

databases$add(
  name = "cyanometdb",
  sourceFiles = list(tsv = "CyanoMetDB_v02_2023.csv"),
  interimFile = "cyanometdb.tsv.gz"
)

databases$add(
  name = "datawarrior",
  sourceFiles = list(tsv = "NaturalProducts.txt"),
  interimFile = "datawarrior.tsv.gz"
)

databases$add(
  name = "dianatdb",
  sourceFiles = list(tsv = "2020_DiaNatDB_336.xlsx"),
  interimFile = "dianatdb.tsv.gz"
)

databases$add(
  name = "dnp",
  sourceFiles = list(tsv = "31_2/full_set.csv"),
  interimFile = "dnp.tsv.gz"
)

databases$add(
  name = "drduke",
  sourceFiles = list(
    tsvCommon = "Duke-Source-CSV/COMMON_NAMES.csv",
    tsvFarmacy = "Duke-Source-CSV/FARMACY_NEW.csv",
    tsvTaxa = "Duke-Source-CSV/FNFTAX.csv",
    tsvReference = "Duke-Source-CSV/REFERENCES.csv"
  ),
  interimFile = "drduke.tsv.gz"
)

databases$add(
  name = "foodb",
  sourceFiles = list(
    tsvCompoundsFlavors = "foodb_2020_04_07_csv/CompoundsFlavor_copy.csv",
    tsvCompounds = "foodb_2020_04_07_csv/Compound_copy.csv",
    tsvContent = "foodb_2020_04_07_csv/Content.csv",
    tsvFlavor = "foodb_2020_04_07_csv/Flavor.csv",
    tsvFood = "foodb_2020_04_07_csv/Food_copy.csv",
    tsvReference = "foodb_2020_04_07_csv/Reference.csv"
  ),
  interimFile = "foodb.tsv.gz"
)

databases$add(
  name = "inflamnat",
  # sourceFiles = list(tsv = "ci8b00560_si_001.xlsx"),
  sourceFiles = list(tsv = "13321_2022_608_MOESM2_ESM.xlsx"),
  interimFile = "inflamnat.tsv.gz"
)

databases$add(
  name = "knapsack",
  sourceFiles = list(tsv = "knapsackScraped.tsv.gz"),
  interimFile = "knapsack.tsv.gz"
)

databases$add(
  name = "custom",
  sourceFiles = NA,
  interimFile = "custom.tsv.gz"
)

databases$add(
  name = "metabolights",
  sourceFiles = list(xmlComplete = "eb-eye_metabolights_complete.xml"),
  interimFile = "metabolights.tsv.gz"
)

databases$add(
  name = "mibig",
  sourceFiles = list(data = "mibig_json_3.1.zip"),
  interimFile = "mibig.tsv.gz"
)

databases$add(
  name = "mitishamba",
  sourceFiles = list(tsv = "mitishambaScraped.tsv.gz"),
  interimFile = "mitishamba.tsv.gz"
)

databases$add(
  name = "napralert",
  sourceFiles = list(
    tsvMatched = "napralert_matched_final_unified.tsv.gz",
    tsvOriginal = "napralert.tsv.gz"
  ),
  interimFile = "napralert.tsv.gz"
)

databases$add(
  name = "npass",
  sourceFiles = list(
    tsvGeneral = "NPASSv2.0_download_naturalProducts_generalInfo.txt",
    tsvSpeciesInfo = "NPASSv2.0_download_naturalProducts_speciesInfo.txt",
    tsvSpeciesPair = "NPASSv2.0_download_naturalProducts_species_pair.txt",
    tsvStructure = "NPASSv2.0_download_naturalProducts_structureInfo.txt"
  ),
  interimFile = "npass.tsv.gz"
)

databases$add(
  name = "npatlas",
  sourceFiles = list(tsv = "NPAtlas_download.tsv"),
  interimFile = "npatlas.tsv.gz"
)

databases$add(
  name = "npcare",
  sourceFiles = list(tsv = "npcare.zip"),
  interimFile = "npcare.tsv.gz"
)

databases$add(
  name = "npedia",
  sourceFiles = list(tsv = "npediaScraped.tsv.gz"),
  interimFile = "npedia.tsv.gz"
)

pathDataExternalDbSourceNubbe <-
  file.path(
    pathDataExternalDbSource,
    "nubbe"
  )

databases$add(
  name = "nubbe",
  sourceFiles = list(tsvPath = file.path(
    "nubbe_raw",
    list.files(
      path = file.path(
        pathDataExternalDbSourceNubbe,
        "nubbe_raw"
      ),
      pattern = "*.xml",
      full.names = FALSE
    )
  )),
  interimFile = "nubbe.tsv.gz"
)

databases$add(
  name = "pamdb",
  sourceFiles = list(tsv = "PaMet.xlsx"),
  interimFile = "pamdb.tsv.gz"
)

databases$add(
  name = "phenolexplorer",
  sourceFiles = list(
    tsvCompoundsClassification = "compounds-classification.csv",
    tsvCompoundsStructures = "compounds-structures.csv",
    tsvCompounds = "compounds.csv",
    tsvFoodsClassification = "foods-classification.csv",
    tsvFoods = "foods.csv",
    tsvMetabolitesStructures = "metabolites-structures.csv",
    tsvMetabolites = "metabolites.csv",
    tsvPublications = "publications.csv",
    tsvComposition = "composition-data.xlsx"
  ),
  interimFile = "phenolexplorer.tsv.gz"
)

databases$add(
  name = "phytohub",
  sourceFiles = list(tsv = "phytohubScraped.tsv.gz"),
  interimFile = "phytohub.tsv.gz"
)

databases$add(
  name = "procardb",
  sourceFiles = list(tsv = "procardbScraped.tsv.gz"),
  interimFile = "procardb.tsv.gz"
)

pathDataExternalDbSourceRespect <-
  file.path(
    pathDataExternalDbSource,
    "respect"
  )

databases$add(
  name = "respect",
  sourceFiles = list(tsvPath = file.path(
    "respect",
    list.files(
      path = file.path(
        pathDataExternalDbSourceRespect,
        "respect"
      ),
      pattern = "*.txt",
      full.names = FALSE
    )
  )),
  interimFile = "respect.tsv.gz"
)

databases$add(
  name = "sancdb",
  sourceFiles = list(tsv = "sancdbScraped.tsv.gz"),
  interimFile = "sancdb.tsv.gz"
)

databases$add(
  name = "streptomedb",
  sourceFiles = list(
    sdf = "streptomedb.sdf",
    tsv = "streptomedb.tsv.gz"
  ),
  interimFile = "streptomedb.tsv.gz"
)

pathDataExternalDbSourceSwmd <-
  file.path(
    pathDataExternalDbSource,
    "swmd"
  )

pathDataExternalDbSourceSwmdDirectory <-
  file.path(
    pathDataExternalDbSourceSwmd,
    "Mol"
  )

databases$add(
  name = "swmd",
  sourceFiles = list(tsv = "swmdScraped.tsv.gz"),
  interimFile = "swmd.tsv.gz"
)

databases$add(
  name = "tmdb",
  sourceFiles = list(tsv = "tmdbScraped.tsv.gz"),
  interimFile = "tmdb.tsv.gz"
)

databases$add(
  name = "tmmc",
  sourceFiles = list(tsv = "compound.xlsx"),
  interimFile = "tmmc.tsv.gz"
)

databases$add(
  name = "tppt",
  sourceFiles = list(tsv = "TPPT_database.xlsx"),
  interimFile = "tppt.tsv.gz"
)

pathDataExternalDbSourceUnpd <-
  file.path(
    pathDataExternalDbSource,
    "unpd"
  )

pathDataExternalDbSourceUnpdIntegrated <-
  file.path(
    pathDataExternalDbSourceUnpd,
    "unpdIntegrated.tsv.gz"
  )

databases$add(
  name = "unpd",
  sourceFiles = list(
    tsvJo = "unpd_final.csv.zip",
    tsvPm = "UNPD_DB.csv.zip"
  ),
  interimFile = "unpd.tsv.gz"
)

databases$add(
  name = "wakankensaku",
  sourceFiles = list(tsv = "wakankensakuScraped.tsv.gz"),
  interimFile = "wakankensaku.tsv.gz"
)

databases$add(
  name = "wikidata",
  sourceFiles = list(tsv = "wikidata.tsv.gz"),
  interimFile = "wikidata.tsv.gz"
)

#### translation source
pathDataExternalTranslationSource <-
  file.path(
    pathDataExternal,
    "translationSource"
  )

##### pubmed
pathDataExternalTranslationSourcePubmed <-
  file.path(
    pathDataExternalTranslationSource,
    "pubmed"
  )

###### file
pathDataExternalTranslationSourcePubmedFile <-
  file.path(
    pathDataExternalTranslationSourcePubmed,
    "PMC-ids.csv.gz"
  )

##### common
pathDataExternalTranslationSourceCommon <-
  file.path(
    pathDataExternalTranslationSource,
    "common"
  )

##### COMMENT: Discrepancy here, don't know if has to be changed #####
###### phenolexplorer
pathDataExternalDbSourcePhenolexplorer <-
  file.path(
    pathDataExternalDbSource,
    "phenolexplorer"
  )

pathDataExternalTranslationSourceCommonPhenolexplorer <-
  file.path(
    pathDataExternalDbSourcePhenolexplorer,
    "foods.csv"
  )

###### foodb
pathDataExternalDbSourceFoodb <-
  file.path(
    pathDataExternalDbSource,
    "foodb"
  )
###### foodb
pathDataExternalTranslationSourceCommonFoodb <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Food_copy.csv"
  )

###### drduke
pathDataExternalDbSourceDrduke <-
  file.path(
    pathDataExternalDbSource,
    "drduke"
  )

###### drduke
####### common
pathDataExternalTranslationSourceCommonDrdukeCommon <-
  file.path(
    pathDataExternalDbSourceDrduke,
    "Duke-Source-CSV/COMMON_NAMES.csv"
  )

####### scientific
pathDataExternalTranslationSourceCommonDrdukeScientific <-
  file.path(
    pathDataExternalDbSourceDrduke,
    "Duke-Source-CSV/FNFTAX.csv"
  )

###### gbif
pathDataExternalTranslationSourceCommonGbif <-
  file.path(
    pathDataExternalTranslationSourceCommon,
    "backbone-current.zip"
  )

##### tcm
pathDataExternalTranslationSourceTcm <-
  file.path(
    pathDataExternalTranslationSource,
    "tcm"
  )

##### COMMENT: Discrepancy here, don't know if has to be changed #####
###### TCMID
pathDataExternalTranslationSourceTcmTcmid <-
  file.path(
    pathDataExternalTranslationSourceTcm,
    "tcmid/data/herb-TCMID.v2.01.txt"
  )

###### Chinese Medicine Board of Australia
pathDataExternalTranslationSourceTcmCmba <-
  file.path(
    pathDataExternalTranslationSourceTcm,
    "Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX"
  )

#### dir
if (mode != "custom") {
  pathDataInterimDbDir <-
    Sys.glob(file.path(paste0(pathDataInterimDb, "/*.tsv.gz")))
} else {
  pathDataInterimDbDir <-
    file.path(pathDataInterimDb, "custom.tsv.gz")
}

#### dictionaries
pathDataInterimDictionaries <- file.path(
  pathDataInterim,
  paste("dictionaries", mode, sep = "_")
)

pathDataInterimDictionariesFix <-
  file.path(
    pathDataExternal,
    "dictionarySource"
  )

##### common
pathDataInterimDictionariesCommon <-
  file.path(
    pathDataInterimDictionariesFix,
    "common"
  )

###### deny
pathDataInterimDictionariesCommonDenyDic <-
  file.path(
    pathDataInterimDictionariesCommon,
    "deny.tsv"
  )

###### manual subtraction
pathDataInterimDictionariesCommonManualSubtraction <-
  file.path(
    pathDataInterimDictionariesCommon,
    "manualSubtraction.tsv"
  )

###### names
pathDataInterimDictionariesCommonNames <-
  file.path(
    pathDataInterimDictionariesCommon,
    "names.tsv.gz"
  )

##### latin
pathDataInterimDictionariesLatin <-
  file.path(
    pathDataInterimDictionariesFix,
    "latin"
  )

###### genitive
pathDataInterimDictionariesLatinGenitive <-
  file.path(
    pathDataInterimDictionariesLatin,
    "genitive"
  )

####### I
pathDataInterimDictionariesLatinGenitiveI <-
  file.path(
    pathDataInterimDictionariesLatinGenitive,
    "i.tsv"
  )

####### Is
pathDataInterimDictionariesLatinGenitiveIs <-
  file.path(
    pathDataInterimDictionariesLatinGenitive,
    "is.tsv"
  )

###### plant parts
pathDataInterimDictionariesLatinPlantParts <-
  file.path(
    pathDataInterimDictionariesLatin,
    "plantParts.tsv"
  )

##### taxa
pathDataInterimDictionariesTaxa <-
  file.path(
    pathDataInterimDictionariesFix,
    "taxa"
  )

###### family
pathDataInterimDictionariesTaxaFamily <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "family.tsv"
  )

###### kingdom
pathDataInterimDictionariesTaxaKingdom <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "kingdom.tsv"
  )

##### organism
pathDataInterimDictionariesOrganism <-
  file.path(
    pathDataInterimDictionaries,
    "organism"
  )

pathDataInterimDictionariesOrganismDictionary <-
  file.path(
    pathDataInterimDictionariesOrganism,
    "organism_dictionary.tsv.gz"
  )

pathDataInterimDictionariesOrganismDictionaryOTL <-
  file.path(
    pathDataInterimDictionariesOrganism,
    "otl.sqlite"
  )

pathDataInterimDictionariesOrganismMetadata <-
  file.path(pathDataInterimDictionariesOrganism, "organism_metadata.tsv.gz")

##### structure
pathDataInterimDictionariesStructure <-
  file.path(
    pathDataInterimDictionaries,
    "structure"
  )

pathDataInterimDictionariesStructureDictionary <-
  file.path(
    pathDataInterimDictionariesStructure,
    "structure_dictionary.tsv.gz"
  )

pathDataInterimDictionariesStructureAntiDictionary <-
  file.path(
    pathDataInterimDictionariesStructure,
    "structure_antiDictionary.tsv.gz"
  )

pathDataInterimDictionariesStructureMetadata <-
  file.path(pathDataInterimDictionariesStructure, "structure_metadata.tsv.gz")

pathDataInterimDictionariesStructureDictionaryChebi <-
  file.path(
    pathDataInterimDictionariesStructure,
    "chebi"
  )

pathDataInterimDictionariesStructureDictionaryChebiFile <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryChebi,
    "chebi.tsv.gz"
  )

pathDataInterimDictionariesStructureDictionaryClassyfire <-
  file.path(
    pathDataInterimDictionariesStructure,
    "classyfire"
  )

pathDataInterimDictionariesStructureDictionaryClassyfireAlternativeParent <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryClassyfire,
    "alternative_parents.tsv.gz"
  )

pathDataInterimDictionariesStructureDictionaryClassyfireDB <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryClassyfire,
    "classyfire.sqlite"
  )

pathDataInterimDictionariesStructureDictionaryClassyfireDirectParent <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryClassyfire,
    "direct_parent.tsv.gz"
  )

pathDataInterimDictionariesStructureDictionaryNpclassifier <-
  file.path(
    pathDataInterimDictionariesStructure,
    "npclassifier"
  )

pathDataInterimDictionariesStructureDictionaryNpclassifierFile <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryNpclassifier,
    "smiles_np_classified.tsv.gz"
  )

##### reference
pathDataInterimDictionariesReference <-
  file.path(
    pathDataInterimDictionaries,
    "reference"
  )

pathDataInterimDictionariesReferenceDictionary <-
  file.path(
    pathDataInterimDictionariesReference,
    "reference_dictionary.tsv.gz"
  )

pathDataInterimDictionariesReferenceOrganismDictionary <-
  file.path(
    pathDataInterimDictionariesReference,
    "reference_dictionaryOrganism.tsv.gz"
  )

pathDataInterimDictionariesReferenceMetadata <-
  file.path(pathDataInterimDictionariesReference, "reference_metadata.tsv.gz")


##### COMMENT: Discrepancy here, don't know if has to be changed #####
####### manual subtraction
pathDataInterimDictionariesTaxaManualSubtraction <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "manualSubtraction.tsv"
  )

####### phylum
pathDataInterimDictionariesTaxaPhylum <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "phylum.tsv"
  )

####### problematic
pathDataInterimDictionariesTaxaProblematic <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "problematic.tsv.zip"
  )

####### wrong verified
pathDataInterimDictionariesTaxaWrongVerified <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "wrongVerified.tsv"
  )

####### wrong homonyms
pathDataInterimDictionariesTaxaWrongHomonyms <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "wrongHomonyms.tsv"
  )

####### ranks
pathDataInterimDictionariesTaxaRanks <-
  file.path(
    pathDataInterimDictionariesTaxa,
    "ranks.tsv"
  )

###### tcm
pathDataInterimDictionariesTcm <-
  file.path(
    pathDataInterimDictionariesFix,
    "tcm"
  )

####### manual subtraction
pathDataInterimDictionariesTcmManualSubtraction <-
  file.path(
    pathDataInterimDictionariesTcm,
    "manualSubtraction.tsv"
  )

####### names
pathDataInterimDictionariesTcmNames <-
  file.path(
    pathDataInterimDictionariesTcm,
    "names.tsv.gz"
  )

#### tables
pathDataInterimTables <- file.path(
  pathDataInterim,
  paste("tables", mode, sep = "_")
)

#### tables
pathDataProcessedTables <- file.path(
  pathDataProcessed,
  paste("tables", mode, sep = "_")
)

#### figures
pathDataProcessedFigures <- file.path(
  pathDataProcessed,
  paste("figures", mode, sep = "_")
)

##### html
pathDataProcessedFiguresHtml <-
  file.path(
    pathDataProcessedFigures,
    "html"
  )

# original fields
pathDataInterimTablesOriginal <-
  file.path(pathDataInterimTables, "0_original")

## structure
pathDataInterimTablesOriginalStructure <-
  file.path(pathDataInterimTablesOriginal, "structure")

### InChI
pathDataInterimTablesOriginalStructureInchi <-
  file.path(pathDataInterimTablesOriginalStructure, "inchi.tsv.gz")

### SMILES
pathDataInterimTablesOriginalStructureSmiles <-
  file.path(pathDataInterimTablesOriginalStructure, "smiles.tsv.gz")

### nominal
pathDataInterimTablesOriginalStructureNominal <-
  file.path(pathDataInterimTablesOriginalStructure, "nominal.tsv.gz")

pathDataInterimTablesOriginalStructureFull <-
  file.path(pathDataInterimTablesOriginalStructure, "full.tsv.gz")

## distinct organism
pathDataInterimTablesOriginalOrganism <-
  file.path(pathDataInterimTablesOriginal, "organism")

### file
pathDataInterimTablesOriginalOrganismFile <-
  file.path(pathDataInterimTablesOriginalOrganism, "original.tsv.gz")

pathDataInterimTablesOriginalOrganismFull <-
  file.path(pathDataInterimTablesOriginalOrganism, "full.tsv.gz")

## ref
pathDataInterimTablesOriginalReference <-
  file.path(pathDataInterimTablesOriginal, "reference")

### DOI
pathDataInterimTablesOriginalReferenceDoi <-
  file.path(pathDataInterimTablesOriginalReference, "doi.tsv.gz")

### original
pathDataInterimTablesOriginalReferenceOriginalFolder <-
  file.path(pathDataInterimTablesOriginalReference, "original")

### pubmed
pathDataInterimTablesOriginalReferencePubmed <-
  file.path(pathDataInterimTablesOriginalReference, "pubmed.tsv.gz")

### title
pathDataInterimTablesOriginalReferencePublishingDetailsFolder <-
  file.path(pathDataInterimTablesOriginalReference, "publishingDetails")

### title
pathDataInterimTablesOriginalReferenceSplitFolder <-
  file.path(pathDataInterimTablesOriginalReference, "split")

### title
pathDataInterimTablesOriginalReferenceTitleFolder <-
  file.path(pathDataInterimTablesOriginalReference, "title")

### full
pathDataInterimTablesOriginalReferenceFull <-
  file.path(pathDataInterimTablesOriginalReference, "full.tsv.gz")

## table
pathDataInterimTablesOriginalTable <-
  file.path(pathDataInterimTablesOriginal, "table.tsv.gz")

# translated fields
pathDataInterimTablesTranslated <-
  file.path(pathDataInterimTables, "1_translated")

## organism
pathDataInterimTablesTranslatedOrganism <-
  file.path(pathDataInterimTablesTranslated, "organism")

### file
##### maybe not useful #####
pathDataInterimTablesTranslatedOrganismFile <-
  file.path(pathDataInterimTablesTranslatedOrganism, "organism.tsv.gz")

## ref
pathDataInterimTablesTranslatedReference <-
  file.path(pathDataInterimTablesTranslated, "reference")

pathDataInterimTablesTranslatedReferenceDoi <-
  file.path(pathDataInterimTablesTranslatedReference, "doi.tsv.gz")

pathDataInterimTablesTranslatedReferenceOriginalFolder <-
  file.path(pathDataInterimTablesTranslatedReference, "original")

pathDataInterimTablesTranslatedReferenceOriginal <-
  file.path(pathDataInterimTablesTranslatedReference, "original.tsv.gz")

pathDataInterimTablesTranslatedReferencePubmed <-
  file.path(pathDataInterimTablesTranslatedReference, "pubmed.tsv.gz")

pathDataInterimTablesTranslatedReferencePublishingDetailsFolder <-
  file.path(
    pathDataInterimTablesTranslatedReference,
    "publishingDetails"
  )

pathDataInterimTablesTranslatedReferencePublishingDetails <-
  file.path(
    pathDataInterimTablesTranslatedReference,
    "publishingDetails.tsv.gz"
  )

pathDataInterimTablesTranslatedReferenceTitleFolder <-
  file.path(pathDataInterimTablesTranslatedReference, "title")

pathDataInterimTablesTranslatedReferenceTitle <-
  file.path(pathDataInterimTablesTranslatedReference, "title.tsv.gz")

pathDataInterimTablesTranslatedReferenceSplitFolder <-
  file.path(pathDataInterimTablesTranslatedReference, "split")

pathDataInterimTablesTranslatedReferenceSplit <-
  file.path(pathDataInterimTablesTranslatedReference, "split.tsv.gz")

pathDataInterimTablesTranslatedReferenceFile <-
  file.path(
    pathDataInterimTablesTranslatedReference,
    "integrated.tsv.gz"
  )

## structure
pathDataInterimTablesTranslatedStructure <-
  file.path(pathDataInterimTablesTranslated, "structure")

### smiles
pathDataInterimTablesTranslatedStructureInchi <-
  file.path(pathDataInterimTablesTranslatedStructure, "inchi.tsv.gz")

### nominal
pathDataInterimTablesTranslatedStructureNominal <-
  file.path(pathDataInterimTablesTranslatedStructure, "nominal.tsv.gz")

#### nominal_opsin
pathDataInterimTablesTranslatedStructureNominal_opsin <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "nominal_opsin.tsv.gz"
  )

#### nominal_pubchem
pathDataInterimTablesTranslatedStructureNominal_pubchem <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "nominal_pubchem.tsv.gz"
  )

### prepared_1
pathDataInterimTablesTranslatedStructurePrepared_1 <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "prepared_1.txt"
  )

### opsin
pathDataInterimTablesTranslatedStructureOpsin <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "opsin.txt"
  )

### both
pathDataInterimTablesTranslatedStructureFinal <-
  file.path(pathDataInterimTablesTranslatedStructure, "final.tsv.gz")

### unique
pathDataInterimTablesTranslatedStructureUnique <-
  file.path(pathDataInterimTablesTranslatedStructure, "unique.tsv.gz")

## table
pathDataInterimTablesTranslatedTable <-
  file.path(pathDataInterimTablesTranslated, "table.tsv.gz")

# processed fields
pathDataInterimTablesProcessed <-
  file.path(pathDataInterimTables, "2_processed")

## organism
pathDataInterimTablesProcessedOrganism <-
  file.path(pathDataInterimTablesProcessed, "organism")

### original
pathDataInterimTablesProcessedOrganismOriginal <-
  file.path(pathDataInterimTablesProcessedOrganism, "original")

pathDataInterimTablesProcessedOrganismOriginalTable <-
  file.path(pathDataInterimTablesProcessedOrganism, "original.tsv.gz")

pathDataInterimTablesProcessedOrganismOriginalUniqueTable <-
  file.path(
    pathDataInterimTablesProcessedOrganism,
    "originalUnique.tsv.gz"
  )

pathDataInterimTablesProcessedOrganismOriginalVerifiedTable <-
  file.path(
    pathDataInterimTablesProcessedOrganism,
    "originalVerified.tsv.gz"
  )

### translated
pathDataInterimTablesProcessedOrganismTranslated <-
  file.path(pathDataInterimTablesProcessedOrganism, "translated")

pathDataInterimTablesProcessedOrganismTranslatedInterim <-
  file.path(pathDataInterimTablesProcessedOrganism, "interim.tsv.gz")

pathDataInterimTablesProcessedOrganismTranslatedTable <-
  file.path(pathDataInterimTablesProcessedOrganism, "translated.tsv.gz")

pathDataInterimTablesProcessedOrganismVerifyTable <-
  file.path(pathDataInterimTablesProcessedOrganism, "verify.tsv.gz")

pathDataInterimTablesProcessedOrganismVerifiedOriginalTable <-
  file.path(
    pathDataInterimTablesProcessedOrganism,
    "original_verified.json"
  )

pathDataInterimTablesProcessedOrganismVerifiedTable <-
  file.path(pathDataInterimTablesProcessedOrganism, "verified.json")

### final processed organisms
pathDataInterimTablesProcessedOrganismFinal <-
  file.path(pathDataInterimTablesProcessedOrganism, "processed.tsv.gz")

pathDataInterimTablesProcessedOrganismRealDiff <-
  file.path(
    pathDataInterimTablesProcessedOrganism,
    "organismsDifferentSpecies.tsv.gz"
  )

### structure
pathDataInterimTablesProcessedStructure <-
  file.path(pathDataInterimTablesProcessed, "structure")

pathDataInterimTablesProcessedStructureFile <-
  file.path(pathDataInterimTablesProcessedStructure, "processed.tsv.gz")

pathDataInterimTablesProcessedStructureStereoCounted <-
  file.path(pathDataInterimTablesProcessedStructure, "counted.tsv.gz")

pathDataInterimTablesProcessedStructureSmiles <-
  file.path(pathDataInterimTablesProcessedStructure, "smiles.tsv.gz")

pathDataInterimTablesProcessedStructureSmiles_1 <-
  file.path(pathDataInterimTablesProcessedStructure, "smiles_1.tsv")

pathDataInterimTablesProcessedStructureSmiles_2 <-
  file.path(pathDataInterimTablesProcessedStructure, "smiles_2.tsv")

pathDataInterimTablesProcessedStructureSmiles_3 <-
  file.path(pathDataInterimTablesProcessedStructure, "smiles_3.tsv")

pathDataInterimTablesProcessedStructureSmiles_4 <-
  file.path(pathDataInterimTablesProcessedStructure, "smiles_4.tsv")

pathDataInterimTablesProcessedStructureNamed <-
  file.path(pathDataInterimTablesProcessedStructure, "named.tsv.gz")

## ref
pathDataInterimTablesProcessedReference <-
  file.path(pathDataInterimTablesProcessed, "reference")

pathDataInterimTablesProcessedReferenceFile <-
  file.path(pathDataInterimTablesProcessedReference, "processed.tsv.gz")

# curated fields
pathDataInterimTablesCurated <-
  file.path(pathDataInterimTables, "3_curated")

### final processed table
pathDataInterimTablesCuratedTable <-
  file.path(pathDataInterimTablesCurated, "table.tsv.gz")

### final processed table
pathDataInterimTablesCuratedTableMaximal <-
  file.path(pathDataInterimTablesCurated, "tableMaximal.tsv.gz")

# analyzed fields
pathDataInterimTablesAnalyzed <-
  file.path(pathDataInterimTables, "4_analyzed")

## triplets
### open
pathDataInterimTablesAnalyzedOpenDbTriplets <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "openDbTriplets.tsv.gz"
  )

### inhouse
pathDataInterimTablesAnalyzedInhouseDbTriplets <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "inhouseDbTriplets.tsv.gz"
  )

### DNP
pathDataInterimTablesAnalyzedClosedDbTriplets <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "closed.tsv.gz"
  )

## structures by kingdom
pathDataInterimTablesAnalyzedStructuresByKingdom <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "structuresByKingdom.tsv"
  )

## unique structures by species
pathDataInterimTablesAnalyzedUniqueStructuresBySpecies <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "uniqueStructuresBySpecies.tsv"
  )

## widespread structures
pathDataInterimTablesAnalyzedWidespreadStructures <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "widespreadStructures.tsv"
  )

## mismatched genera
pathDataInterimTablesAnalyzedMismatchedGenera <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "mismatchedGenera.tsv"
  )

## redundancy table
pathDataInterimTablesAnalyzedRedundancyTable <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "redundancyTable.tsv"
  )

## sample ONPDB triplets (all)
pathDataInterimTablesAnalyzedSampleAllONPDB <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "sampleAllONPDB.tsv"
  )

## validated structure-organism pairs
pathDataInterimTablesAnalyzedPlatinum <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "validated_referenced_structure_organism_pairs.tsv.gz"
  )

pathDataInterimTablesAnalyzedPlatinumNew <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "validated_referenced_structure_organism_pairs_new.tsv.gz"
  )

## sample knapsack triplets
pathDataInterimTablesAnalyzedSampleKnapsack <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "sampleKnapsack.tsv"
  )

pathDataInterimTimestamps <-
  file.path(
    pathDataInterim,
    "timestamps.txt"
  )

## dirty for the moment
pathOriginalGnfinderScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_full.sh",
  "custom" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_custom.sh",
  "min" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_test.sh"
)

pathTranslatedGnfinderScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_full.sh",
  "custom" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_custom.sh",
  "min" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_test.sh"
)

pathOriginalGnverifierScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_full.sh",
  "custom" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_custom.sh",
  "min" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_test.sh"
)

pathGnverifierScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/gnverifierLauncher_full.sh",
  "custom" = "2_curating/2_editing/organism/shell/gnverifierLauncher_custom.sh",
  "min" = "2_curating/2_editing/organism/shell/gnverifierLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/gnverifierLauncher_test.sh"
)

pathTimestampsScript <- "3_analyzing/timestamps.sh"

pathTests <- file.path("../tests")

pathTestsFile <-
  file.path(pathTests, switch(mode,
    "full" = "tests_min.tsv",
    "min" = "tests_min.tsv",
    "test" = "tests.tsv"
  ))

pathTestsDicCommonFile <-
  file.path(pathTests, "tests_dic_common_min.tsv.gz")

pathTestsDicTcmFile <-
  file.path(pathTests, "tests_dic_tcm_min.tsv.gz")

pathTestsExpectations <- file.path(pathTests, "expectations")

pathTestsOrganisms <- file.path(
  pathTestsExpectations,
  "organisms.tsv"
)

pathTestsStructures <- file.path(
  pathTestsExpectations,
  "structures.tsv"
)

pathTestsReferences <- file.path(
  pathTestsExpectations,
  "references.tsv"
)

pathTestsPlatinum <- file.path(
  pathTestsExpectations,
  "validated.tsv"
)

path_accepted_fields <- file.path(
  "1_gathering",
  "accepted_fields.tsv"
)

pathLastWdExport <- "../data/interim/db/wikidata.tsv.gz"

pathLastFrozen <- "230106_frozen_metadata.csv.gz"

pathLastFrozenClosed <- "230106_closed_metadata.csv.gz"

pathDataInterimTablesAnalyzedGarbage <-
  file.path(
    pathDataInterimTablesAnalyzed,
    "recycle.tsv.gz"
  )

pathLastTreeBio <- "tree_bio.json"

pathLastTreeChemo <- "tree_chemo.json"
