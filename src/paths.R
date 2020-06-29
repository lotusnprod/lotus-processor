#######################################################
######################   Paths   ######################
#######################################################

source("functions/database.R")

# root
## data
pathData <- Sys.getenv("DATA_PATH",
                       unset = "../data")


### external
pathDataExternal <-
  file.path(pathData,
            "external")

#### db source
pathDataExternalDbSource <-
  file.path(pathDataExternal,
            "dbSource")

### interim
# pathDataInterim <-
#   file.path(pathData,
#             "interim")

pathDataInterim <-
  file.path(pathData,
            "interim_min")

#### db
pathDataInterimDb <-
  file.path(pathDataInterim,
            "db")


databases <-
  Databases$new(pathDbSource = pathDataExternalDbSource,
                pathDbInterim = pathDataInterimDb)

databases$add(
  name = "afrotryp",
  sourceFiles = list(tsv = "afrotryp.tsv.zip"),
  interimFile = "afrotryp.tsv.zip"
)

databases$add(
  name = "alkamid",
  sourceFiles = list(tsv = "alkamidScraped.tsv.zip",
                     tsvRef = "alkamidRefScraped.tsv.zip"),
  interimFile = "alkamid.tsv.zip"
)

databases$add(
  name = "biofacquim",
  sourceFiles = list(sdf = "BIOFACQUIM_V2.sdf",
                     tsv = "biofacquim.tsv.zip"),
  interimFile = "biofacquim.tsv.zip"
)

databases$add(
  name = "biophytmol",
  sourceFiles = list(tsv = "biophytmolScraped.tsv.zip"),
  interimFile = "biophytmol.tsv.zip"
)

databases$add(
  name = "carotenoiddb",
  sourceFiles = list(tsv = "carotenoiddbScraped.tsv.zip",
                     tsvInchi = "Carotenoids_InChI_InChIKey.tsv"),
  interimFile = "carotenoiddb.tsv.zip"
)

databases$add(
  name = "cmaup",
  sourceFiles = list(
    tsvIngredients = "CMAUPv1.0_download_Ingredients_All.txt",
    tsvPlants = "CMAUPv1.0_download_Plants.txt",
    tsvAssociations = "CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt"
  ),
  interimFile = "cmaup.tsv.zip"
)

databases$add(
  name = "coconut",
  sourceFiles = list(sdf = "COCONUT_DB.sdf",
                     tsv = "coconutConverted.tsv.zip"),
  interimFile = "coconut.tsv.zip"
)

databases$add(
  name = "cyanometdb",
  sourceFiles = list(tsv = "media-1.csv"),
  interimFile = "cyanometdb.tsv.zip"
)

databases$add(
  name = "dnp",
  sourceFiles = list(tsv = "28_2/full_set.csv"),
  interimFile = "dnp.tsv.zip"
)

###### COMMENT THERE ARE SCRIPTS IN DNP DATA FOLDER #######
###### WE HAVE TO MOVE THEM WHEN UPDATING #######

databases$add(
  name = "drduke",
  sourceFiles = list(
    tsvCommon = "Duke-Source-CSV/COMMON_NAMES.csv",
    tsvFarmacy = "Duke-Source-CSV/FARMACY_NEW.csv",
    tsvTaxa = "Duke-Source-CSV/FNFTAX.csv",
    tsvReference = "Duke-Source-CSV/REFERENCES.csv"
  ),
  interimFile = "drduke.tsv.zip"
)

# COMMENT not sure about how clean those lines are
##### etcm
pathDataExternalDbSourceEtcm <-
  file.path(pathDataExternalDbSource,
            "etcm")

databases$add(
  name = "etcm",
  sourceFiles = list(tsvPath = file.path(
    "data",
    list.files(
      path = file.path(pathDataExternalDbSourceEtcm,
                       "data"),
      pattern = "*.csv",
      full.names = FALSE
    )
  )),
  interimFile = "etcm.tsv.zip"
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
  interimFile = "foodb.tsv.zip"
)

databases$add(
  name = "inflamnat",
  sourceFiles = list(tsv = "ci8b00560_si_001.xlsx"),
  interimFile = "inflamnat.tsv.zip"
)

databases$add(
  name = "knapsack",
  sourceFiles = list(tsv = "knapsackScraped.tsv.zip"),
  interimFile = "knapsack.tsv.zip"
)

databases$add(
  name = "metabolights",
  sourceFiles = list(
    xmlComplete = "eb-eye_metabolights_complete.xml",
    xmlStudies = "eb-eye_metabolights_studies.xml",
    tsvPrecleaned = "metabolightsPrecleaned.tsv.zip",
    tsvStudies = "metabolightsStudiesScraped.tsv.zip"
  ),
  interimFile = "metabolights.tsv.zip"
)

###### COMMENT Not clean #######
pathDataExternalDbSourceMetabolights <-
  file.path(pathDataExternalDbSource,
            "metabolights")

###### studies scraped directory
pathDataExternalDbSourceMetabolightsStudiesScrapedDir <-
  file.path(pathDataExternalDbSourceMetabolights,
            "studiesScraped")

##### mibig
pathDataExternalDbSourceMibig <-
  file.path(pathDataExternalDbSource,
            "mibig")

databases$add(
  name = "mibig",
  sourceFiles = list(tsvPath = file.path(
    "mibig_json_2.0",
    list.files(
      path = file.path(pathDataExternalDbSourceMibig,
                       "mibig_json_2.0"),
      pattern = "*.json",
      full.names = FALSE
    )
  )),
  interimFile = "mibig.tsv.zip"
)

databases$add(
  name = "mitishamba",
  sourceFiles = list(tsv = "mitishambaScraped.tsv.zip"),
  interimFile = "mitishamba.tsv.zip"
)

databases$add(
  name = "nanpdb",
  sourceFiles = list(tsv = "nanpdbScraped.tsv.zip"),
  interimFile = "nanpdb.tsv.zip"
)

databases$add(
  name = "npass",
  sourceFiles = list(
    tsvGeneral = "NPASSv1.0_download_naturalProducts_generalInfo.txt",
    tsvProperties = "NPASSv1.0_download_naturalProducts_properties.txt",
    tsvSpeciesInfo = "NPASSv1.0_download_naturalProducts_speciesInfo.txt",
    tsvSpeciesPair = "NPASSv1.0_download_naturalProducts_species_pair.txt"
  ),
  interimFile = "npass.tsv.zip"
)

databases$add(
  name = "npatlas",
  sourceFiles = list(tsv = "np_atlas_2019_12.tsv"),
  interimFile = "npatlas.tsv.zip"
)

databases$add(
  name = "npcare",
  sourceFiles = list(tsv = "npcare.zip"),
  interimFile = "npcare.tsv.zip"
)

databases$add(
  name = "npedia",
  sourceFiles = list(tsv = "npediaScraped.tsv.zip"),
  interimFile = "npedia.tsv.zip"
)

databases$add(
  name = "pamdb",
  sourceFiles = list(tsv = "PaMet.xlsx"),
  interimFile = "pamdb.tsv.zip"
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
  interimFile = "phenolexplorer.tsv.zip"
)

databases$add(
  name = "phytohub",
  sourceFiles = list(tsv = "phytohubScraped.tsv.zip"),
  interimFile = "phytohub.tsv.zip"
)

# here we should allow a list of interim files, so for the moment I let it as it is
# databases$add(
#   name = "plantcyc",
#   sourceFiles = list(tsv = "phytohubScraped.tsv.zip"),
#   interimFile = "plantcyc.tsv.zip"
# )

##### plantcyc
pathDataExternalDbSourcePlantcyc <-
  file.path(pathDataExternalDbSource,
            "plantcyc")

###### directory
pathDataExternalDbSourcePlantcycDir <-
  file.path(pathDataExternalDbSourcePlantcyc,
            "0_data")

###### original
pathDataExternalDbSourcePlantcycOriginal <-
  list.files(path = pathDataExternalDbSourcePlantcyc,
             pattern = "*.tsv.zip",
             full.names = TRUE)

databases$add(
  name = "procardb",
  sourceFiles = list(tsv = "procardbScraped.tsv.zip"),
  interimFile = "procardb.tsv.zip"
)

##### respect
pathDataExternalDbSourceRespect <-
  file.path(pathDataExternalDbSource,
            "respect")

databases$add(
  name = "respect",
  sourceFiles = list(tsvPath = file.path(
    "respect",
    list.files(
      path = file.path(pathDataExternalDbSourceRespect,
                       "respect"),
      pattern = "*.txt",
      full.names = FALSE
    )
  )),
  interimFile = "respect.tsv.zip"
)

databases$add(
  name = "sancdb",
  sourceFiles = list(tsv = "sancdbScraped.tsv.zip"),
  interimFile = "sancdb.tsv.zip"
)

databases$add(
  name = "streptomedb",
  sourceFiles = list(sdf = "streptomedb.sdf",
                     tsv = "streptomedb.tsv.zip"),
  interimFile = "streptomedb.tsv.zip"
)

##### swmd
pathDataExternalDbSourceSwmd <-
  file.path(pathDataExternalDbSource,
            "swmd")

##### directory
pathDataExternalDbSourceSwmdDirectory <-
  file.path(pathDataExternalDbSourceSwmd,
            "Mol")

databases$add(
  name = "swmd",
  sourceFiles = list(tsv = "swmdScraped.tsv.zip"),
  interimFile = "swmd.tsv.zip"
)

##### symmap
pathDataExternalDbSourceSymmap <-
  file.path(pathDataExternalDbSource,
            "symmap")

databases$add(
  name = "symmap",
  sourceFiles = list(tsvPath = file.path(
    "data",
    list.files(
      path = file.path(pathDataExternalDbSourceSymmap,
                       "data"),
      pattern = "*.csv",
      full.names = FALSE
    )
  )),
  interimFile = "symmap.tsv.zip"
)

###### bio
pathDataExternalDbSourceSymmapBio <-
  file.path(pathDataExternalDbSourceSymmap,
            "SymMap v1.0, SMHB file.xlsx")

###### chemo
pathDataExternalDbSourceSymmapChemo <-
  file.path(pathDataExternalDbSourceSymmap,
            "SymMap v1.0, SMIT file.xlsx")

databases$add(
  name = "tmdb",
  sourceFiles = list(tsv = "tmdbScraped.tsv.zip"),
  interimFile = "tmdb.tsv.zip"
)

databases$add(
  name = "tmmc",
  sourceFiles = list(tsv = "compound.xlsx"),
  interimFile = "tmmc.tsv.zip"
)

databases$add(
  name = "tppt",
  sourceFiles = list(tsv = "TPPT_database.xlsx"),
  interimFile = "tppt.tsv.zip"
)

databases$add(
  name = "triforc",
  sourceFiles = list(tsv1 = "triforcOriginal.tsv",
                     tsv2 = "triforcBis.tsv"),
  interimFile = "triforc.tsv.zip"
)

##### triforc
pathDataExternalDbSourceTriforc <-
  file.path(pathDataExternalDbSource,
            "triforc")

##### to get
pathDataExternalDbSourceTriforcToGet <-
  file.path(pathDataExternalDbSourceTriforc,
            "triforcToGet.tsv")

databases$add(
  name = "triforc",
  sourceFiles = list(tsv1 = "triforcOriginal.tsv",
                     tsv2 = "triforcBis.tsv"),
  interimFile = "triforc.tsv.zip"
)

##### unpd
pathDataExternalDbSourceUnpd <-
  file.path(pathDataExternalDbSource,
            "unpd")

##### compiled
pathDataExternalDbSourceUnpdIntegrated <-
  file.path(pathDataExternalDbSourceUnpd,
            "unpdIntegrated.tsv.zip")

databases$add(
  name = "unpd",
  sourceFiles = list(tsvJo = "unpd_final.csv.zip",
                     tsvPm = "UNPD_DB.csv.zip"),
  interimFile = "unpd.tsv.zip"
)

#### translation source
pathDataExternalTranslationSource <-
  file.path(pathDataExternal,
            "translationSource")

##### common
pathDataExternalTranslationSourceCommon <-
  file.path(pathDataExternalTranslationSource,
            "common")

##### COMMENT: Discrepancy here, don't know if has to be changed #####

###### phenolexplorer

###### phenolexplorer
pathDataExternalDbSourcePhenolexplorer <-
  file.path(pathDataExternalDbSource,
            "phenolexplorer")

pathDataExternalTranslationSourceCommonPhenolexplorer <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "foods.csv")

###### foodb
pathDataExternalDbSourceFoodb <-
  file.path(pathDataExternalDbSource,
            "foodb")
###### foodb
pathDataExternalTranslationSourceCommonFoodb <-
  file.path(pathDataExternalDbSourceFoodb,
            "foodb_2020_04_07_csv/Food_copy.csv")

###### drduke
pathDataExternalDbSourceDrduke <-
  file.path(pathDataExternalDbSource,
            "drduke")

###### drduke
####### common
pathDataExternalTranslationSourceCommonDrdukeCommon <-
  file.path(pathDataExternalDbSourceDrduke,
            "Duke-Source-CSV/COMMON_NAMES.csv")

####### scientific
pathDataExternalTranslationSourceCommonDrdukeScientific <-
  file.path(pathDataExternalDbSourceDrduke,
            "Duke-Source-CSV/FNFTAX.csv")

###### gbif
####### vernacular
pathDataExternalTranslationSourceCommonGbifVernacular <-
  file.path(
    pathDataExternalTranslationSourceCommon,
    "backbone-current/VernacularName.tsv.zip"
  )

####### scientific
pathDataExternalTranslationSourceCommonGbifScientific <-
  file.path(pathDataExternalTranslationSourceCommon,
            "backbone-current/Taxon.tsv.zip")

##### tcm
pathDataExternalTranslationSourceTcm <-
  file.path(pathDataExternalTranslationSource,
            "tcm")

##### COMMENT: Discrepancy here, don't know if has to be changed #####

###### TCMID
pathDataExternalTranslationSourceTcmTcmid <-
  file.path(pathDataExternalTranslationSourceTcm,
            "TCMID/data/herb-TCMID.v2.01.txt")

###### Chinese Medicine Board of Australia
pathDataExternalTranslationSourceTcmCmba <-
  file.path(
    pathDataExternalTranslationSourceTcm,
    "Chinese-Medicine-Board---List---Nomenclature-list-of-commonly-used-Chinese-herbal-medicines.XLSX"
  )

#### dir
pathDataInterimDbDir <-
  Sys.glob(file.path(paste(pathDataInterimDb,
                           "/*.tsv.zip",
                           sep = "")))


##### biofacquim
pathDataInterimDbBiofacquim <-
  file.path(pathDataInterimDb,
            "biofacquim.tsv.zip")

##### biophytmol
pathDataInterimDbBiophytmol <-
  file.path(pathDataInterimDb,
            "biophytmol.tsv.zip")

##### carotenoiddb
pathDataInterimDbCarotenoiddb <-
  file.path(pathDataInterimDb,
            "carotenoiddb.tsv.zip")

##### cmaup
pathDataInterimDbCmaup <-
  file.path(pathDataInterimDb,
            "cmaup.tsv.zip")

##### coconut
pathDataInterimDbCoconut <-
  file.path(pathDataInterimDb,
            "coconut.tsv.zip")

##### cyanometdb
pathDataInterimDbCyanometdb <-
  file.path(pathDataInterimDb,
            "cyanometdb.tsv.zip")

##### dnp
pathDataInterimDbDnp <-
  file.path(pathDataInterimDb,
            "dnp.tsv.zip")

##### drduke
pathDataInterimDbDrduke <-
  file.path(pathDataInterimDb,
            "drduke.tsv.zip")

##### etcm
pathDataInterimDbEtcm <-
  file.path(pathDataInterimDb,
            "etcm.tsv.zip")

##### foodb
pathDataInterimDbFoodb <-
  file.path(pathDataInterimDb,
            "foodb.tsv.zip")

##### inflamnat
pathDataInterimDbInflamnat <-
  file.path(pathDataInterimDb,
            "inflamnat.tsv.zip")

##### knapsack
pathDataInterimDbKnapsack <-
  file.path(pathDataInterimDb,
            "knapsack.tsv.zip")

##### metabolights
pathDataInterimDbMetabolights <-
  file.path(pathDataInterimDb,
            "metabolights.tsv.zip")

##### mibig
pathDataInterimDbMibig <-
  file.path(pathDataInterimDb,
            "mibig.tsv.zip")

##### mitishamba
pathDataInterimDbMitishamba <-
  file.path(pathDataInterimDb,
            "mitishamba.tsv.zip")

##### nanpdb
pathDataInterimDbNanpdb <-
  file.path(pathDataInterimDb,
            "nanpdb.tsv.zip")

##### npass
pathDataInterimDbNpass <-
  file.path(pathDataInterimDb,
            "npass.tsv.zip")

##### npatlas
pathDataInterimDbNpatlas <-
  file.path(pathDataInterimDb,
            "npatlas.tsv.zip")

##### npcare
pathDataInterimDbNpcare <-
  file.path(pathDataInterimDb,
            "npcare.tsv.zip")

##### npedia
pathDataInterimDbNpedia <-
  file.path(pathDataInterimDb,
            "npedia.tsv.zip")

##### pamdb
pathDataInterimDbPamdb <-
  file.path(pathDataInterimDb,
            "pamdb.tsv.zip")

##### phenolexplorer
pathDataInterimDbPhenolexplorer <-
  file.path(pathDataInterimDb,
            "phenolexplorer.tsv.zip")

##### phytohub
pathDataInterimDbPhytohub <-
  file.path(pathDataInterimDb,
            "phytohub.tsv.zip")

##### plantcyc
pathDataInterimDbPlantcyc <-
  file.path(pathDataInterimDb,
            "plantcyc.tsv.zip")

##### procardb
pathDataInterimDbProcardb <-
  file.path(pathDataInterimDb,
            "procardb.tsv.zip")

##### respect
pathDataInterimDbRespect <-
  file.path(pathDataInterimDb,
            "respect.tsv.zip")

##### sancdb
pathDataInterimDbSancdb <-
  file.path(pathDataInterimDb,
            "sancdb.tsv.zip")

##### streptomedb
pathDataInterimDbStreptomedb <-
  file.path(pathDataInterimDb,
            "streptomedb.tsv.zip")

##### swmd
pathDataInterimDbSwmd <-
  file.path(pathDataInterimDb,
            "swmd.tsv.zip")

##### symmap
pathDataInterimDbSymmap <-
  file.path(pathDataInterimDb,
            "symmap.tsv.zip")

##### tmdb
pathDataInterimDbTmdb <-
  file.path(pathDataInterimDb,
            "tmdb.tsv.zip")

##### tmmc
pathDataInterimDbTmmc <-
  file.path(pathDataInterimDb,
            "tmmc.tsv.zip")

##### tppt
pathDataInterimDbTmmc <-
  file.path(pathDataInterimDb,
            "tppt.tsv.zip")

##### triforc
pathDataInterimDbTriforc <-
  file.path(pathDataInterimDb,
            "triforc.tsv.zip")

##### unpd
pathDataInterimDbUnpd <-
  file.path(pathDataInterimDb,
            "unpd.tsv.zip")

#### dictionaries
pathDataInterimDictionaries <-
  file.path(pathDataInterim,
            "dictionaries")

###### common
pathDataInterimDictionariesCommon <-
  file.path(pathDataInterimDictionaries,
            "common")

####### black
pathDataInterimDictionariesCommonBlackDic <-
  file.path(pathDataInterimDictionariesCommon,
            "black.tsv")

####### manual subtraction
pathDataInterimDictionariesCommonManualSubtraction <-
  file.path(pathDataInterimDictionariesCommon,
            "manualSubtraction.tsv")

####### names
pathDataInterimDictionariesCommonNames <-
  file.path(pathDataInterimDictionariesCommon,
            "names.tsv.zip")

###### latin
pathDataInterimDictionariesLatin <-
  file.path(pathDataInterimDictionaries,
            "latin")

####### genitive
pathDataInterimDictionariesLatinGenitive <-
  file.path(pathDataInterimDictionariesLatin,
            "genitive")

######## I
pathDataInterimDictionariesLatinGenitiveI <-
  file.path(pathDataInterimDictionariesLatinGenitive,
            "i.tsv")

######## Is
pathDataInterimDictionariesLatinGenitiveIs <-
  file.path(pathDataInterimDictionariesLatinGenitive,
            "is.tsv")

####### plant parts
pathDataInterimDictionariesLatinPlantParts <-
  file.path(pathDataInterimDictionariesLatin,
            "plantParts.tsv")

###### taxa
pathDataInterimDictionariesTaxa <-
  file.path(pathDataInterimDictionaries,
            "taxa")

####### family
pathDataInterimDictionariesTaxaFamily <-
  file.path(pathDataInterimDictionariesTaxa,
            "family.tsv")

####### kingdom
pathDataInterimDictionariesTaxaKingdom <-
  file.path(pathDataInterimDictionariesTaxa,
            "kingdom.tsv")

##### COMMENT: Discrepancy here, don't know if has to be changed #####

####### manual subtraction
pathDataInterimDictionariesTaxaManualSubtraction <-
  file.path(pathDataInterimDictionariesTaxa,
            "manualSubtraction.tsv")

####### phylum
pathDataInterimDictionariesTaxaPhylum <-
  file.path(pathDataInterimDictionariesTaxa,
            "phylum.tsv")

####### problematic
pathDataInterimDictionariesTaxaProblematic <-
  file.path(pathDataInterimDictionariesTaxa,
            "problematic.tsv.zip")

####### problematic
pathDataInterimDictionariesTaxaRanks <-
  file.path(pathDataInterimDictionariesTaxa,
            "ranks.tsv")

###### tcm
pathDataInterimDictionariesTcm <-
  file.path(pathDataInterimDictionaries,
            "tcm")

####### manual subtraction
pathDataInterimDictionariesTcmManualSubtraction <-
  file.path(pathDataInterimDictionariesTcm,
            "manualSubtraction.tsv")

####### names
pathDataInterimDictionariesTcmNames <-
  file.path(pathDataInterimDictionariesTcm,
            "names.tsv.zip")

##### COMMENT ##### The tables will have to slowly move to lists in dictionaries folder in flat format.s

#### tables
pathDataInterimTables <-
  file.path(pathDataInterim,
            "tables")

### processed
pathDataProcessed <-
  file.path(pathData,
            "processed")

#### figures
pathDataProcessedFigures <-
  file.path(pathDataProcessed,
            "figures")

##### html
pathDataProcessedFiguresHtml <-
  file.path(pathDataProcessedFigures,
            "html")




# original fields
## structure
### InChI
pathOriginalStructureInchi <-
  file.path(pathDataInterimTables, "0_original/originalStructureInchi.tsv.zip")

### SMILES
pathOriginalStructureSmiles <- 
  file.path(pathDataInterimTables, "0_original/originalStructureSmiles.tsv.zip")

### nominal
pathOriginalStructureNominal <-
  "../data/interim/tables/0_original/originalStructureNominal.tsv.zip"

## organism
pathOriginalOrganism <-
  "../data/interim/tables/0_original/originalOrganism.tsv.zip"

## distinct organism
pathOriginalOrganismDistinct <-
  "../data/interim/tables/0_original/gnfinder/"

## ref
pathOriginalRef <-
  "../data/interim/tables/0_original/originalReference.tsv.zip"

## table
pathOriginalTable <-
  "../data/interim/tables/0_original/originalTable.tsv.zip"

# translated fields
## organism
pathTranslatedOrganism <-
  "../data/interim/tables/1_translated/translatedOrganism.tsv.zip"

## ref
pathTranslatedReference <-
  "../data/interim/tables/1_translated/translatedReference.tsv.zip"

## structure
### smiles
pathTranslatedStructureSmiles <-
  "../data/interim/tables/1_translated/translatedStructureSmiles.tsv.zip"

### nominal
pathTranslatedStructureNominal <-
  "../data/interim/tables/1_translated/translatedStructureNominal.tsv.zip"

### nominal_dvc
pathTranslatedStructureNominal_dvc <-
  "../data/dvc_pipeline_outputs/tables/1_translated/translatedStructureNominal.tsv.zip"


## distinct organism
pathTranslatedOrganismDistinct <-
  "../data/interim/tables/1_translated/gnfinder/"

## distinct structure
pathTranslatedStructureDistinct <-
  "../data/interim/tables/1_translated/rdkit/translatedStructureRdkit.tsv.zip"

## table
pathTranslatedTable <-
  "../data/interim/tables/1_translated/translatedTable.tsv.zip"

# cleaned fields
## ref
pathCleanedReference <-
  "../data/interim/tables/2_cleaned/cleanedReference.tsv.zip"

## organism
### gnfinder
#### original
##### json dir
pathCleanedOrganismOriginalDirJson <-
  "../data/interim/tables/2_cleaned/gnfinder/original/json/"

##### tsv converted dir
pathCleanedOrganismOriginalDirTsv <-
  "../data/interim/tables/2_cleaned/gnfinder/original/tsv/"

#### translated
##### json dir
pathCleanedOrganismTranslatedDirJson <-
  "../data/interim/tables/2_cleaned/gnfinder/translated/json/"

##### tsv converted dir
pathCleanedOrganismTranslatedDirTsv <-
  "../data/interim/tables/2_cleaned/gnfinder/translated/tsv/"

### final cleaned organisms
pathCleanedOrganism <-
  "../data/interim/tables/2_cleaned/cleanedOrganism.tsv.zip"

pathCuratedOrganism <-
  "../data/interim/tables/3_curated/curatedOrganism.tsv.zip"

pathCuratedOrganismRealDiff <-
  "../data/interim/tables/curatedOrganismsDifferentSpecies.tsv.zip"

## dirty for the moment
pathOriginalGnfinderScript <-
  "2_curating/2_editing/bio/subscripts/shell/originalGnfinderLauncher.sh"

pathTranslatedGnfinderScript <-
  "2_curating/2_editing/bio/subscripts/shell/translatedGnfinderLauncher.sh"
