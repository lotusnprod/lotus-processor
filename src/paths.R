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
pathDataInterim <-
  file.path(pathData,
            "interim")

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
  sourceFiles = list(tsv = "apps_database_csv_BIOFACQUIM.csv"),
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
  sourceFiles = list(sdf = "COCONUT.sdf.zip",
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
  sourceFiles = list(tsvCommon = "Duke-Source-CSV/COMMON_NAMES.csv",
                     tsvFarmacy = "Duke-Source-CSV/FARMACY_NEW.csv",
                     tsvTaxa = "Duke-Source-CSV/FNFTAX.csv",
                     tsvReference = "Duke-Source-CSV/REFERENCES.csv"),
  interimFile = "drduke.tsv.zip"
)

# COMMENT not sure about those lines
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

##### metabolights
pathDataExternalDbSourceMetabolights <-
  file.path(pathDataExternalDbSource,
            "metabolights")

###### COMMENT Not clean #######

###### complete
pathDataExternalDbSourceMetabolightsComplete <-
  file.path(pathDataExternalDbSourceMetabolights,
            "eb-eye_metabolights_complete.xml")

###### studies
pathDataExternalDbSourceMetabolightsStudies <-
  file.path(pathDataExternalDbSourceMetabolights,
            "eb-eye_metabolights_studies.xml")

###### studies scraped directory
pathDataExternalDbSourceMetabolightsStudiesScrapedDir <-
  file.path(pathDataExternalDbSourceMetabolights,
            "studiesScraped")

###### precleaned
pathDataExternalDbSourceMetabolightsPrecleaned <-
  file.path(pathDataExternalDbSourceMetabolights,
            "metabolightsPrecleaned.tsv.zip")

###### studies scraped
pathDataExternalDbSourceMetabolightsStudiesScraped <-
  file.path(pathDataExternalDbSourceMetabolights,
            "metabolightsStudiesScraped.tsv.zip")

##### mibig
pathDataExternalDbSourceMibig <-
  file.path(pathDataExternalDbSource,
            "mibig")

###### original
pathDataExternalDbSourceMibigOriginal <-
  list.files(
    path = file.path(pathDataExternalDbSourceMibig,
                     "mibig_json_2.0/"),
    pattern = "*.json",
    full.names = TRUE
  )

##### mitishamba
pathDataExternalDbSourceMitishamba <-
  file.path(pathDataExternalDbSource,
            "mitishamba")

###### original
pathDataExternalDbSourceMitishambaOriginal <-
  file.path(pathDataExternalDbSourceMitishamba,
            "mitishambaScraped.tsv.zip")

##### nanpdb
pathDataExternalDbSourceNanpdb <-
  file.path(pathDataExternalDbSource,
            "nanpdb")

###### original
pathDataExternalDbSourceNanpdbOriginal <-
  file.path(pathDataExternalDbSourceNanpdb,
            "nanpdbScraped.tsv.zip")

##### npass
pathDataExternalDbSourceNpass <-
  file.path(pathDataExternalDbSource,
            "npass")

###### general info
pathDataExternalDbSourceNpassGeneralInfo <-
  file.path(
    pathDataExternalDbSourceNpass,
    "NPASSv1.0_download_naturalProducts_generalInfo.txt"
  )

###### properties
pathDataExternalDbSourceNpassProperties <-
  file.path(
    pathDataExternalDbSourceNpass,
    "NPASSv1.0_download_naturalProducts_properties.txt"
  )

###### species info
pathDataExternalDbSourceNpassSpeciesInfo <-
  file.path(
    pathDataExternalDbSourceNpass,
    "NPASSv1.0_download_naturalProducts_speciesInfo.txt"
  )

###### species pair
pathDataExternalDbSourceNpassSpeciesPair <-
  file.path(
    pathDataExternalDbSourceNpass,
    "NPASSv1.0_download_naturalProducts_species_pair.txt"
  )

##### npatlas
pathDataExternalDbSourceNpatlas <-
  file.path(pathDataExternalDbSource,
            "npatlas")

###### original
pathDataExternalDbSourceNpatlasOriginal <-
  file.path(pathDataExternalDbSourceNpatlas,
            "np_atlas_2019_12.tsv")

##### npcare
pathDataExternalDbSourceNpcare <-
  file.path(pathDataExternalDbSource,
            "npcare")

###### original
pathDataExternalDbSourceNpcareOriginal <-
  file.path(pathDataExternalDbSourceNpcare,
            "npcare.zip")

##### npedia
pathDataExternalDbSourceNpedia <-
  file.path(pathDataExternalDbSource,
            "npedia")

###### original
pathDataExternalDbSourceNpediaOriginal <-
  file.path(pathDataExternalDbSourceNpedia,
            "npediaScraped.tsv.zip")

##### pamdb
pathDataExternalDbSourcePamdb <-
  file.path(pathDataExternalDbSource,
            "pamdb")

###### original
pathDataExternalDbSourcePamdbOriginal <-
  file.path(pathDataExternalDbSourcePamdb,
            "PaMet.xlsx")

##### phenolexplorer
pathDataExternalDbSourcePhenolexplorer <-
  file.path(pathDataExternalDbSource,
            "phenolexplorer")

###### compounds classification
pathDataExternalDbSourcePhenolexplorerCompoundsClassification <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "compounds-classification.csv")

###### compounds structure
pathDataExternalDbSourcePhenolexplorerCompoundsStructures <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "compounds-structures.csv")

###### compounds
pathDataExternalDbSourcePhenolexplorerCompounds <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "compounds.csv")

###### foods classification
pathDataExternalDbSourcePhenolexplorerFoodsClassification <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "foods-classification.csv")

###### foods
pathDataExternalDbSourcePhenolexplorerFoods <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "foods.csv")

###### metabolites structure
pathDataExternalDbSourcePhenolexplorerMetabolitesStructures <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "metabolites-structures.csv")

###### metabolites
pathDataExternalDbSourcePhenolexplorerMetabolites <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "metabolites.csv")

###### publications
pathDataExternalDbSourcePhenolexplorerPublications <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "publications.csv")

###### composition
pathDataExternalDbSourcePhenolexplorerComposition <-
  file.path(pathDataExternalDbSourcePhenolexplorer,
            "composition-data.xlsx")

##### phytohub
pathDataExternalDbSourcePhytohub <-
  file.path(pathDataExternalDbSource,
            "phytohub")

###### original
pathDataExternalDbSourcePhytohubOriginal <-
  file.path(pathDataExternalDbSourcePhytohub,
            "phytohubScraped.tsv.zip")

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

##### procardb
pathDataExternalDbSourceProcardb <-
  file.path(pathDataExternalDbSource,
            "procardb")

###### original
pathDataExternalDbSourceProcardbOriginal <-
  file.path(pathDataExternalDbSourceProcardb,
            "procardbScraped.tsv.zip")

##### respect
pathDataExternalDbSourceRespect <-
  file.path(pathDataExternalDbSource,
            "respect")

###### directory
pathDataExternalDbSourceRespectDir <-
  file.path(pathDataExternalDbSourceRespect,
            "respect")

##### sancdb
pathDataExternalDbSourceSancdb <-
  file.path(pathDataExternalDbSource,
            "sancdb")

###### original
pathDataExternalDbSourceSancdbOriginal <-
  file.path(pathDataExternalDbSourceSancdb,
            "sancdbScraped.tsv.zip")

##### streptomedb
pathDataExternalDbSourceStreptomedb <-
  file.path(pathDataExternalDbSource,
            "streptomedb")

##### original
pathDataExternalDbSourceStreptomedbOriginal <-
  file.path(pathDataExternalDbSourceStreptomedb,
            "streptomedb.sdf")

##### compiled
pathDataExternalDbSourceStreptomedbCompiled <-
  file.path(pathDataExternalDbSourceStreptomedb,
            "streptomedb.tsv.zip")

##### swmd
pathDataExternalDbSourceSwmd <-
  file.path(pathDataExternalDbSource,
            "swmd")

##### directory
pathDataExternalDbSourceSwmdDirectory <-
  file.path(pathDataExternalDbSourceSwmd,
            "Mol")

##### original
pathDataExternalDbSourceSwmdOriginal <-
  file.path(pathDataExternalDbSourceSwmd,
            "swmdScraped.tsv.zip")

##### symmap
pathDataExternalDbSourceSymmap <-
  file.path(pathDataExternalDbSource,
            "symmap")

###### original
pathDataExternalDbSourceSymmapOriginal <-
  list.files(
    path = file.path(pathDataExternalDbSourceSymmap,
                     "data"),
    pattern = "*.csv",
    full.names = TRUE
  )

###### bio
pathDataExternalDbSourceSymmapBio <-
  file.path(pathDataExternalDbSourceSymmap,
            "SymMap v1.0, SMHB file.xlsx")

###### chemo
pathDataExternalDbSourceSymmapChemo <-
  file.path(pathDataExternalDbSourceSymmap,
            "SymMap v1.0, SMIT file.xlsx")

##### tmdb
pathDataExternalDbSourceTmdb <-
  file.path(pathDataExternalDbSource,
            "tmdb")

##### original
pathDataExternalDbSourceTmdbOriginal <-
  file.path(pathDataExternalDbSourceTmdb,
            "tmdbScraped.tsv.zip")

##### tmmc
pathDataExternalDbSourceTmmc <-
  file.path(pathDataExternalDbSource,
            "tmmc")

##### original
pathDataExternalDbSourceTmmcOriginal <-
  file.path(pathDataExternalDbSourceTmmc,
            "compound.xlsx")

##### tppt
pathDataExternalDbSourceTppt <-
  file.path(pathDataExternalDbSource,
            "tppt")

##### original
pathDataExternalDbSourceTpptOriginal <-
  file.path(pathDataExternalDbSourceTppt,
            "TPPT_database.xlsx")

##### triforc
pathDataExternalDbSourceTriforc <-
  file.path(pathDataExternalDbSource,
            "triforc")

##### original
pathDataExternalDbSourceTriforcOriginal <-
  file.path(pathDataExternalDbSourceTriforc,
            "triforcOriginal.tsv")

##### to get
pathDataExternalDbSourceTriforcToGet <-
  file.path(pathDataExternalDbSourceTriforc,
            "triforcToGet.tsv")

##### bis
pathDataExternalDbSourceTriforBis <-
  file.path(pathDataExternalDbSourceTriforc,
            "triforcBis.tsv")

##### unpd
pathDataExternalDbSourceUnpd <-
  file.path(pathDataExternalDbSource,
            "unpd")

##### original Jo
pathDataExternalDbSourceUnpdOriginal_1 <-
  file.path(pathDataExternalDbSourceUnpd,
            "unpd_final.csv.zip")

##### original PM
pathDataExternalDbSourceUnpdOriginal_2 <-
  file.path(pathDataExternalDbSourceUnpd,
            "UNPD_DB.csv.zip")

##### compiled
pathDataExternalDbSourceUnpdCompiled <-
  file.path(pathDataExternalDbSourceUnpd,
            "unpdCompiled.tsv.zip")

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
  "../data/interim/tables/0_original/originalStructureInchi.tsv.zip"

### SMILES
pathOriginalStructureSmiles <-
  "../data/interim/tables/0_original/originalStructureSmiles.tsv.zip"

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
