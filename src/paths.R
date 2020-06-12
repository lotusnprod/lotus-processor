#######################################################
######################   Paths   ######################
#######################################################

# root
## data
pathData <- "../data"

### external
pathDataExternal <-
  file.path(pathData,
            "external")

#### db source
pathDataExternalDbSource <-
  file.path(pathDataExternal,
            "dbSource")

##### afrotryp
pathDataExternalDbSourceAfrotryp <-
  file.path(pathDataExternalDbSource,
            "afrotryp")

###### original
pathDataExternalDbSourceAfrotrypOriginal <-
  file.path(pathDataExternalDbSourceAfrotryp,
            "afrotryp.tsv.zip")

##### alkamid
pathDataExternalDbSourceAlkamid <-
  file.path(pathDataExternalDbSource,
            "alkamid")

###### original
pathDataExternalDbSourceAlkamidOriginal <-
  file.path(pathDataExternalDbSourceAlkamid,
            "alkamidScraped.tsv.zip")

###### ref
pathDataExternalDbSourceAlkamidRef <-
  file.path(pathDataExternalDbSourceAlkamid,
            "alkamidRefScraped.tsv.zip")

##### biofacquim
pathDataExternalDbSourceBiofacquim <-
  file.path(pathDataExternalDbSource,
            "biofacquim")

###### original
pathDataExternalDbSourceBiofacquimOriginal <-
  file.path(
    pathDataExternalDbSourceBiofacquim,
    "apps_database_csv_BIOFACQUIM.csv"
  )

##### biophytmol
pathDataExternalDbSourceBiophytmol <-
  file.path(pathDataExternalDbSource,
            "biophytmol")

###### original
pathDataExternalDbSourceBiophytmolOriginal <-
  file.path(pathDataExternalDbSourceBiophytmol,
            "biophytmolScraped.tsv.zip")

##### carotenoiddb
pathDataExternalDbSourceCarotenoiddb <-
  file.path(pathDataExternalDbSource,
            "carotenoiddb")

###### original
pathDataExternalDbSourceCarotenoiddbOriginal <-
  file.path(
    pathDataExternalDbSourceCarotenoiddb,
    "carotenoiddbScraped.tsv.zip"
  )

###### inchikey
pathDataExternalDbSourceCarotenoiddbInchikey <-
  file.path(
    pathDataExternalDbSourceCarotenoiddb,
    "Carotenoids_InChI_InChIKey.tsv"
  )

##### cmaup
pathDataExternalDbSourceCmaup <-
  file.path(pathDataExternalDbSource,
            "cmaup")

###### ingredients
pathDataExternalDbSourceCmaupIngredients <-
  file.path(
    pathDataExternalDbSourceCmaup,
    "CMAUPv1.0_download_Ingredients_All.txt"
  )

###### plants
pathDataExternalDbSourceCmaupPlants <-
  file.path(pathDataExternalDbSourceCmaup,
            "CMAUPv1.0_download_Plants.txt")

###### associations
pathDataExternalDbSourceCmaupAssociations <-
  file.path(
    pathDataExternalDbSourceCmaup,
    "CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt"
  )

##### coconut
pathDataExternalDbSourceCoconut <-
  file.path(pathDataExternalDbSource,
            "coconut")

###### original
pathDataExternalDbSourceCoconutOriginal <-
  file.path(pathDataExternalDbSourceCoconut,
            "coconutCompiled.tsv.zip")

##### cyanometdb
pathDataExternalDbSourceCyanometdb <-
  file.path(pathDataExternalDbSource,
            "cyanometdb")

###### original
pathDataExternalDbSourceCyanometdbOriginal <-
  file.path(pathDataExternalDbSourceCyanometdb,
            "media-1.tsv.zip")

##### dnp
pathDataExternalDbSourceDnp <-
  file.path(pathDataExternalDbSource,
            "dnp")

###### COMMENT THERE ARE SCRIPTS IN DNP DATA FOLDER #######
###### WE HAVE TO MOVE THEM WHEN UPDATING #######

###### original
pathDataExternalDbSourceDnpOriginal <-
  file.path(pathDataExternalDbSourceDnp,
            "28_2/full_set.csv")

##### drduke
pathDataExternalDbSourceDrduke <-
  file.path(pathDataExternalDbSource,
            "drduke")

###### common names
pathDataExternalDbSourceDrdukeCommonNames <-
  file.path(
    pathDataExternalDbSourceDrduke,
    "Duke-Source-CSV/COMMON_NAMES.csv"
  )

###### farmacy
pathDataExternalDbSourceDrdukeFarmacy <-
  file.path(
    pathDataExternalDbSourceDrduke,
    "Duke-Source-CSV/FARMACY_NEW.csv"
  )

###### taxa
pathDataExternalDbSourceDrdukeTaxa <-
  file.path(pathDataExternalDbSourceDrduke,
            "Duke-Source-CSV/FNFTAX.csv")

###### references
pathDataExternalDbSourceDrdukeReferences <-
  file.path(
    pathDataExternalDbSourceDrduke,
    "Duke-Source-CSV/REFERENCES.csv"
  )

##### etcm
pathDataExternalDbSourceEtcm <-
  file.path(pathDataExternalDbSource,
            "etcm")

###### original
pathDataExternalDbSourceEtcmOriginal <-
  list.files(
    path = file.path(pathDataExternalDbSourceEtcm,
                     "data/"),
    pattern = "*.csv",
    full.names = TRUE
  )

##### foodb
pathDataExternalDbSourceFoodb <-
  file.path(pathDataExternalDbSource,
            "foodb")

###### compounds flavors
pathDataExternalDbSourceFoodbCompoundsFlavors <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/CompoundsFlavor_copy.csv"
  )

###### compounds
pathDataExternalDbSourceFoodbCompounds <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Compound_copy.csv"
  )

###### content
pathDataExternalDbSourceFoodbContent <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Content.csv"
  )

###### flavor
pathDataExternalDbSourceFoodbFlavor <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Flavor.csv"
  )

###### food
pathDataExternalDbSourceFoodbFood <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Food_copy.csv"
  )

###### reference
pathDataExternalDbSourceFoodbReference <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Reference.csv"
  )

##### inflamnat
pathDataExternalDbSourceInflamnat <-
  file.path(pathDataExternalDbSource,
            "inflamnat")

###### original
pathDataExternalDbSourceInflamnatOriginal <-
  file.path(pathDataExternalDbSourceInflamnat,
            "ci8b00560_si_001.xlsx")

##### knapsack
pathDataExternalDbSourceKnapsack <-
  file.path(pathDataExternalDbSource,
            "knapsack")

###### original
pathDataExternalDbSourceKnapsackOriginal <-
  file.path(pathDataExternalDbSourceKnapsack,
            "knapsackScraped.tsv.zip")

##### metabolights
pathDataExternalDbSourceMetabolights <-
  file.path(pathDataExternalDbSource,
            "metabolights")

###### COMMENT Not clean #######

###### complete
pathDataExternalDbSourceMetabolightsComplete <-
  file.path(
    pathDataExternalDbSourceMetabolights,
    "eb-eye_metabolights_complete.xml"
  )

###### studies
pathDataExternalDbSourceMetabolightsStudies <-
  file.path(
    pathDataExternalDbSourceMetabolights,
    "eb-eye_metabolights_studies.xml"
  )

###### studies scraped directory
pathDataExternalDbSourceMetabolightsStudiesScrapedDir <-
  file.path(pathDataExternalDbSourceMetabolights,
            "studiesScraped")

###### precleaned
pathDataExternalDbSourceMetabolightsPrecleaned <-
  file.path(
    pathDataExternalDbSourceMetabolights,
    "metabolightsPrecleaned.tsv.zip"
  )

###### studies scraped
pathDataExternalDbSourceMetabolightsStudiesScraped <-
  file.path(
    pathDataExternalDbSourceMetabolights,
    "metabolightsStudiesScraped.tsv.zip"
  )

##### mibig
pathDataExternalDbSourceMibig <-
  file.path(pathDataExternalDbSource,
            "mibig")

###### original
pathDataExternalDbSourceMibigOriginal <-
  list.files(
    path = file.path(
      pathDataExternalDbSourceMibig,
      "mibig_json_2.0/"
    ),
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

##### npedia
pathDataExternalDbSourceNpedia <-
  file.path(pathDataExternalDbSource,
            "npedia")

##### pamdb
pathDataExternalDbSourcePamdb <-
  file.path(pathDataExternalDbSource,
            "pamdb")

##### phenolexplorer
pathDataExternalDbSourcePhenolexplorer <-
  file.path(pathDataExternalDbSource,
            "phenolexplorer")

##### phytohub
pathDataExternalDbSourcePhytohub <-
  file.path(pathDataExternalDbSource,
            "phytohub")

##### plantcyc
pathDataExternalDbSourcePlantcyc <-
  file.path(pathDataExternalDbSource,
            "plantcyc")

##### procardb
pathDataExternalDbSourceProcardb <-
  file.path(pathDataExternalDbSource,
            "procardb")

##### respect
pathDataExternalDbSourceRespect <-
  file.path(pathDataExternalDbSource,
            "respect")

##### sancdb
pathDataExternalDbSourceSancdb <-
  file.path(pathDataExternalDbSource,
            "sancdb")

##### streptomedb
pathDataExternalDbSourceStreptomedb <-
  file.path(pathDataExternalDbSource,
            "streptomedb")

##### swmd
pathDataExternalDbSourceSwmd <-
  file.path(pathDataExternalDbSource,
            "swmd")

##### symmap
pathDataExternalDbSourceSymmap <-
  file.path(pathDataExternalDbSource,
            "symmap")

##### tmdb
pathDataExternalDbSourceTmdb <-
  file.path(pathDataExternalDbSource,
            "tmdb")

##### tmmc
pathDataExternalDbSourceTmmc <-
  file.path(pathDataExternalDbSource,
            "tmmc")

##### tppt
pathDataExternalDbSourceTppt <-
  file.path(pathDataExternalDbSource,
            "tppt")

##### triforc
pathDataExternalDbSourceTriforc <-
  file.path(pathDataExternalDbSource,
            "triforc")

##### unpd
pathDataExternalDbSourceTriforc <-
  file.path(pathDataExternalDbSource,
            "unpd")

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
pathDataExternalTranslationSourceCommonFoodb <-
  file.path(
    pathDataExternalDbSourceFoodb,
    "foodb_2020_04_07_csv/Food_copy.csv"
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
    "DRDUKE/Duke-Source-CSV/FNFTAX.csv"
  )

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

###### TM-MC
pathDataExternalTranslationSourceTcmTmmc <-
  file.path(pathDataExternalDbSourceTmmc,
            "TMMC/compound.xlsx")

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

### interim
pathDataInterim <-
  file.path(pathData,
            "interim")

#### db
pathDataInterimDb <-
  file.path(pathDataInterim,
            "db")

##### afrotryp
pathDataInterimDbAfrotryp <-
  file.path(pathDataInterimDb,
            "afrotryp.tsv.zip")

##### alkamid
pathDataInterimDbAlkamid <-
  file.path(pathDataInterimDb,
            "alkamid.tsv.zip")

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














# dbInterim
## allDBs
pathAllDbsInterim <-
  Sys.glob(file.path("../data/interim/db/*_std.tsv.zip"))

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

# sanitized fields
## ref
pathSanitizedReference <-
  "../data/interim/tables/2_sanitized/sanitizedReference.tsv.zip"

## organism
### gnfinder json dir
pathSanitizedOrganismDirJson <-
  "../data/interim/tables/2_sanitized/gnfinder/json/"

### tsv converted dir
pathSanitizedOrganismDirTsv <-
  "../data/interim/tables/2_sanitized/gnfinder/tsv/"

### final sanitized organisms
pathSanitizedOrganism <-
  "../data/interim/tables/2_sanitized/sanitizedOrganism.tsv.zip"
