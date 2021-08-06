###############################################################################
##################################   Paths   ##################################
###############################################################################

source("r/database.R")

mode <- Sys.getenv("MODE", unset = "full")

if (exists("mode_test")) {
  mode <- "test"
}

works_locally_only <- TRUE

## for steps where paths need to be adapted locally or related to other repos
molconvertPath <- "~/../../Applications/MarvinSuite/bin/molconvert"
wikidataLotusExporterPath <- "../../wikidataLotusExporter"
wikidataLotusExporterDataPath <-
  file.path(wikidataLotusExporterPath, "data")
wikidataLotusExporterDataOutputPath <-
  file.path(wikidataLotusExporterDataPath, "output")
wikidataLotusExporterDataOutputTaxaPath <-
  file.path(wikidataLotusExporterDataOutputPath, "taxa.tsv")
wikidataLotusExporterDataOutputStructuresPath <-
  file.path(wikidataLotusExporterDataOutputPath, "compounds.tsv")
wikidataLotusExporterDataOutputReferencesPath <-
  file.path(wikidataLotusExporterDataOutputPath, "references.tsv")
wikidataLotusExporterDataOutputTriplesPath <-
  file.path(
    wikidataLotusExporterDataOutputPath,
    "compound_reference_taxon.tsv"
  )

# databases for which we have no right to disseminate the data
forbidden_export <- c("dnp", "foodb")

# root
## bin
pathBin <- Sys.getenv("BIN_PATH",
  unset = "../bin"
)

## opsin
pathBinOpsin <- file.path(
  pathBin,
  "opsin-2.5.0-jar-with-dependencies.jar"
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
pathDataInterimDb <-
  file.path(
    pathDataInterim,
    "db"
  )

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
  sourceFiles = list(tsv = "1-s2.0-S0043135421002153-mmc4.zip"),
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
  sourceFiles = list(tsv = "30_1/full_set.csv"),
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
  sourceFiles = list(tsv = "ci8b00560_si_001.xlsx"),
  interimFile = "inflamnat.tsv.gz"
)

databases$add(
  name = "knapsack",
  sourceFiles = list(tsv = "knapsackScraped.tsv.gz"),
  interimFile = "knapsack.tsv.gz"
)

databases$add(
  name = "metabolights",
  sourceFiles = list(xmlComplete = "eb-eye_metabolights_complete.xml"),
  interimFile = "metabolights.tsv.gz"
)

databases$add(
  name = "mibig",
  sourceFiles = list(data = "data.zip"),
  interimFile = "mibig.tsv.gz"
)

databases$add(
  name = "mitishamba",
  sourceFiles = list(tsv = "mitishambaScraped.tsv.gz"),
  interimFile = "mitishamba.tsv.gz"
)

databases$add(
  name = "nanpdb",
  sourceFiles = list(tsv = "nanpdbScraped.tsv.gz"),
  interimFile = "nanpdb.tsv.gz"
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
    tsvGeneral = "NPASSv1.0_download_naturalProducts_generalInfo.txt",
    tsvProperties = "NPASSv1.0_download_naturalProducts_properties.txt",
    tsvSpeciesInfo = "NPASSv1.0_download_naturalProducts_speciesInfo.txt",
    tsvSpeciesPair = "NPASSv1.0_download_naturalProducts_species_pair.txt"
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
    "Chinese-Medicine-Board---List---Nomenclature-list-of-commonly-used-Chinese-herbal-medicines.XLSX"
  )

#### dir
pathDataInterimDbDir <-
  Sys.glob(file.path(paste0(pathDataInterimDb, "/*.tsv.gz")))

#### dictionaries
pathDataInterimDictionaries <- switch(mode,
  "full" = file.path(
    pathDataInterim,
    "dictionaries"
  ),
  "min" = file.path(
    pathDataInterim,
    "dictionaries_min"
  ),
  "test" = file.path(
    pathDataInterim,
    "dictionaries_test"
  )
)

pathDataInterimDictionariesFix <-
  file.path(
    pathDataInterim,
    "dictionaries"
  )

##### common
pathDataInterimDictionariesCommon <-
  file.path(
    pathDataInterimDictionariesFix,
    "common"
  )

###### black
pathDataInterimDictionariesCommonBlackDic <-
  file.path(
    pathDataInterimDictionariesCommon,
    "black.tsv"
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
    "dictionary.tsv.gz"
  )

pathDataInterimDictionariesOrganismDictionaryOTL <-
  file.path(
    pathDataInterimDictionariesOrganism,
    "otl.sqlite"
  )

pathDataInterimDictionariesOrganismMetadata <-
  file.path(pathDataInterimDictionariesOrganism, "metadata.tsv.gz")

##### structure
pathDataInterimDictionariesStructure <-
  file.path(
    pathDataInterimDictionaries,
    "structure"
  )

pathDataInterimDictionariesStructureDictionary <-
  file.path(
    pathDataInterimDictionariesStructure,
    "dictionary.tsv.gz"
  )

pathDataInterimDictionariesStructureAntiDictionary <-
  file.path(
    pathDataInterimDictionariesStructure,
    "antiDictionary.tsv.gz"
  )

pathDataInterimDictionariesStructureMetadata <-
  file.path(pathDataInterimDictionariesStructure, "metadata.tsv.gz")

pathDataInterimDictionariesStructureDictionaryClassyfire <-
  file.path(
    pathDataInterimDictionariesFix,
    "structure/classyfire"
  )

pathDataInterimDictionariesStructureDictionaryClassyfireFile <-
  file.path(
    pathDataInterimDictionariesStructureDictionaryClassyfire,
    "classy.tsv.gz"
  )

pathDataInterimDictionariesStructureDictionaryNpclassifier <-
  file.path(
    pathDataInterimDictionariesFix,
    "structure/npclassifier"
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
    "dictionary.tsv.gz"
  )

pathDataInterimDictionariesReferenceOrganismDictionary <-
  file.path(
    pathDataInterimDictionariesReference,
    "dictionaryOrganism.tsv.gz"
  )

pathDataInterimDictionariesReferenceMetadata <-
  file.path(pathDataInterimDictionariesReference, "metadata.tsv.gz")


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
pathDataInterimTables <- switch(mode,
  "full" = file.path(
    pathDataInterim,
    "tables"
  ),
  "min" = file.path(
    pathDataInterim,
    "tables_min"
  ),
  "test" = file.path(
    pathDataInterim,
    "tables_test"
  )
)

#### tables
pathDataProcessedTables <- switch(mode,
  "full" = file.path(
    pathDataProcessed,
    "tables"
  ),
  "min" = file.path(
    pathDataProcessed,
    "tables_min"
  ),
  "test" = file.path(
    pathDataProcessed,
    "tables_test"
  )
)

#### figures
pathDataProcessedFigures <- switch(mode,
  "full" = file.path(
    pathDataProcessed,
    "figures"
  ),
  "min" = file.path(
    pathDataProcessed,
    "figures_min"
  ),
  "test" = file.path(
    pathDataProcessed,
    "figures_test"
  )
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

### publishingDetails
pathDataInterimTablesOriginalReferencePublishingDetails <-
  file.path(
    pathDataInterimTablesOriginalReference,
    "publishingDetails.tsv.gz"
  )

### title
pathDataInterimTablesOriginalReferenceTitleFolder <-
  file.path(pathDataInterimTablesOriginalReference, "title")

### split
pathDataInterimTablesOriginalReferenceSplit <-
  file.path(pathDataInterimTablesOriginalReference, "split.tsv.gz")

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

pathDataInterimTablesTranslatedReferencePublishingDetails <-
  file.path(
    pathDataInterimTablesTranslatedReference,
    "publishingDetails.tsv.gz"
  )

pathDataInterimTablesTranslatedReferenceTitleFolder <-
  file.path(pathDataInterimTablesTranslatedReference, "title")

pathDataInterimTablesTranslatedReferenceTitle <-
  file.path(pathDataInterimTablesTranslatedReference, "title.tsv.gz")

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
pathDataInterimTablesTranslatedStructureSmiles <-
  file.path(pathDataInterimTablesTranslatedStructure, "smiles.tsv.gz")

### nominal
pathDataInterimTablesTranslatedStructureNominal <-
  file.path(pathDataInterimTablesTranslatedStructure, "nominal.tsv.gz")

#### nominal_opsin
pathDataInterimTablesTranslatedStructureNominal_opsin <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "nominal_opsin.tsv.gz"
  )

#### nominal_cts
pathDataInterimTablesTranslatedStructureNominal_cts <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "nominal_cts.tsv.gz"
  )

#### nominal_cts_2
pathDataInterimTablesTranslatedStructureNominal_cts_2 <-
  file.path(
    pathDataInterimTablesTranslatedStructure,
    "nominal_cts_2.tsv.gz"
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

# cleaned fields
pathDataInterimTablesCleaned <-
  file.path(pathDataInterimTables, "2_cleaned")

## organism
pathDataInterimTablesCleanedOrganism <-
  file.path(pathDataInterimTablesCleaned, "organism")

### original
pathDataInterimTablesCleanedOrganismOriginal <-
  file.path(pathDataInterimTablesCleanedOrganism, "original")

pathDataInterimTablesCleanedOrganismOriginalTable <-
  file.path(pathDataInterimTablesCleanedOrganism, "original.tsv.gz")

pathDataInterimTablesCleanedOrganismOriginalUniqueTable <-
  file.path(
    pathDataInterimTablesCleanedOrganism,
    "originalUnique.tsv.gz"
  )

pathDataInterimTablesCleanedOrganismOriginalVerifiedTable <-
  file.path(
    pathDataInterimTablesCleanedOrganism,
    "originalVerified.tsv.gz"
  )

### translated
pathDataInterimTablesCleanedOrganismTranslated <-
  file.path(pathDataInterimTablesCleanedOrganism, "translated")

pathDataInterimTablesCleanedOrganismTranslatedInterim <-
  file.path(pathDataInterimTablesCleanedOrganism, "interim.tsv.gz")

pathDataInterimTablesCleanedOrganismTranslatedTable <-
  file.path(pathDataInterimTablesCleanedOrganism, "translated.tsv.gz")

pathDataInterimTablesCleanedOrganismVerifyTable <-
  file.path(pathDataInterimTablesCleanedOrganism, "verify.tsv.gz")

pathDataInterimTablesCleanedOrganismVerifiedOriginalTable <-
  file.path(
    pathDataInterimTablesCleanedOrganism,
    "original_verified.json"
  )

pathDataInterimTablesCleanedOrganismVerifiedTable <-
  file.path(pathDataInterimTablesCleanedOrganism, "verified.json")

### final cleaned organisms
pathDataInterimTablesCleanedOrganismFinal <-
  file.path(pathDataInterimTablesCleanedOrganism, "cleaned.tsv.gz")

pathDataInterimTablesCleanedOrganismRealDiff <-
  file.path(
    pathDataInterimTablesCleanedOrganism,
    "organismsDifferentSpecies.tsv.gz"
  )

### structure
pathDataInterimTablesCleanedStructure <-
  file.path(pathDataInterimTablesCleaned, "structure")

pathDataInterimTablesCleanedStructureFile <-
  file.path(pathDataInterimTablesCleanedStructure, "cleaned.tsv.gz")

pathDataInterimTablesCleanedStructureStereoCounted <-
  file.path(pathDataInterimTablesCleanedStructure, "counted.tsv.gz")

pathDataInterimTablesCleanedStructureSmiles <-
  file.path(pathDataInterimTablesCleanedStructure, "smiles.tsv.gz")

pathDataInterimTablesCleanedStructureSmiles_1 <-
  file.path(pathDataInterimTablesCleanedStructure, "smiles_1.tsv")

pathDataInterimTablesCleanedStructureSmiles_2 <-
  file.path(pathDataInterimTablesCleanedStructure, "smiles_2.tsv")

pathDataInterimTablesCleanedStructureSmiles_3 <-
  file.path(pathDataInterimTablesCleanedStructure, "smiles_3.tsv")

pathDataInterimTablesCleanedStructureSmiles_4 <-
  file.path(pathDataInterimTablesCleanedStructure, "smiles_4.tsv")

pathDataInterimTablesCleanedStructureNamed <-
  file.path(pathDataInterimTablesCleanedStructure, "named.tsv.gz")

## ref
pathDataInterimTablesCleanedReference <-
  file.path(pathDataInterimTablesCleaned, "reference")

pathDataInterimTablesCleanedReferenceFile <-
  file.path(pathDataInterimTablesCleanedReference, "cleaned.tsv.gz")

# curated fields
pathDataInterimTablesCurated <-
  file.path(pathDataInterimTables, "3_curated")

### final cleaned table
pathDataInterimTablesCuratedTable <-
  file.path(pathDataInterimTablesCurated, "table.tsv.gz")

### final cleaned table
pathDataInterimTablesCuratedTableMaximal <-
  file.path(pathDataInterimTablesCurated, "tableMaximal.tsv.gz")

# analysed fields
pathDataInterimTablesAnalysed <-
  file.path(pathDataInterimTables, "4_analysed")

## triplets
### open
pathDataInterimTablesAnalysedOpenDbTriplets <-
  file.path(
    pathDataInterimTablesAnalysed,
    "openDbTriplets.tsv.gz"
  )

### inhouse
pathDataInterimTablesAnalysedInhouseDbTriplets <-
  file.path(
    pathDataInterimTablesAnalysed,
    "inhouseDbTriplets.tsv.gz"
  )

### DNP
pathDataInterimTablesAnalysedDnpDbTriplets <-
  file.path(
    pathDataInterimTablesAnalysed,
    "dnp.tsv.gz"
  )

## structures by kingdom
pathDataInterimTablesAnalysedStructuresByKingdom <-
  file.path(
    pathDataInterimTablesAnalysed,
    "structuresByKingdom.tsv"
  )

## unique structures by species
pathDataInterimTablesAnalysedUniqueStructuresBySpecies <-
  file.path(
    pathDataInterimTablesAnalysed,
    "uniqueStructuresBySpecies.tsv"
  )

## widespread structures
pathDataInterimTablesAnalysedWidespreadStructures <-
  file.path(
    pathDataInterimTablesAnalysed,
    "widespreadStructures.tsv"
  )

## mismatched genera
pathDataInterimTablesAnalysedMismatchedGenera <-
  file.path(
    pathDataInterimTablesAnalysed,
    "mismatchedGenera.tsv"
  )

## redundancy table
pathDataInterimTablesAnalysedRedundancyTable <-
  file.path(
    pathDataInterimTablesAnalysed,
    "redundancyTable.tsv"
  )

## sample ONPDB triplets (all)
pathDataInterimTablesAnalysedSampleAllONPDB <-
  file.path(
    pathDataInterimTablesAnalysed,
    "sampleAllONPDB.tsv"
  )

## sample ONPDB triplets (gold)
pathDataInterimTablesAnalysedGold <-
  file.path(
    pathDataInterimTablesAnalysed,
    "gold.tsv.gz"
  )

pathDataInterimTablesAnalysedPlatinum <-
  file.path(
    pathDataInterimTablesAnalysed,
    "platinum.tsv.gz"
  )

## sample ONPDB triplets (gold)
pathDataInterimTablesAnalysedSampleGoldONPDB <-
  file.path(
    pathDataInterimTablesAnalysed,
    "sampleGoldONPDB.tsv"
  )

## sample knapsack triplets
pathDataInterimTablesAnalysedSampleKnapsack <-
  file.path(
    pathDataInterimTablesAnalysed,
    "sampleKnapsack.tsv"
  )

## dirty for the moment
pathOriginalGnfinderScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_full.sh",
  "min" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/originalGnfinderLauncher_test.sh"
)

pathTranslatedGnfinderScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_full.sh",
  "min" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/translatedGnfinderLauncher_test.sh"
)

pathOriginalGnverifierScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_full.sh",
  "min" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/originalGnverifierLauncher_test.sh"
)

pathGnverifierScript <- switch(mode,
  "full" = "2_curating/2_editing/organism/shell/gnverifierLauncher_full.sh",
  "min" = "2_curating/2_editing/organism/shell/gnverifierLauncher_min.sh",
  "test" = "2_curating/2_editing/organism/shell/gnverifierLauncher_test.sh"
)

pathTests <- file.path("../tests")

pathTestsFile <- file.path(pathTests, "tests.tsv")

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

path_accepted_fields <- file.path(
  "1_gathering",
  "accepted_fields.tsv"
)

pathLastWdExport <- "210505_wikidata_query.tsv"

pathLastFrozen <- "210715_frozen_metadata.csv.gz"

pathLastFrozenDnp <- "210715_dnp_metadata.csv.gz"
