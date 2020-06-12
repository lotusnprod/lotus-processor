#######################################################
######################   Paths   ######################
#######################################################

# translationSource
## common
### phenolexplorer
pathTranslationSourceCommonPhenolexplorer <-
  "../data/external/dbSource/PHENOLEXPLORER/0_initial_files/foods.csv"

### foodb
pathTranslationSourceCommonFoodb <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Food_copy.csv"

### drduke
#### common
pathTranslationSourceCommonDrdukeCommon <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/COMMON_NAMES.csv"

#### scientific
pathTranslationSourceCommonDrdukeScientific <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/FNFTAX.csv"

### gbif
#### vernacular
pathTranslationSourceCommonGbifVernacular <-
  "../data/external/translationSource/common/backbone-current/VernacularName.tsv.zip"

#### scientific
pathTranslationSourceCommonGbifScientific <-
  "../data/external/translationSource/common/backbone-current/Taxon.tsv.zip"

# interim (common names)
## manual subtraction
pathInterimCommonManualSubtraction <-
  "../data/interim/dictionaries/commonManualSubtraction.tsv"

## common names dic
pathInterimCommonNamesDic <-
  "../data/interim/dictionaries/commonNamesDic.tsv.zip"

## tcm
### TM-MC
pathTranslationSourceTcmTmmc <-
  "../data/external/dbSource/TMMC/0_initial_files/compound.xlsx"

### TCMID
pathTranslationSourceTcmTcmid <-
  "../data/external/translationSource/tcm/TCMID/data/herb-TCMID.v2.01.txt"

### Chinese Medicine Board of Australia
pathTranslationSourceTcmCmba <-
  "../data/external/translationSource/tcm/Chinese-Medicine-Board---List---Nomenclature-list-of-commonly-used-Chinese-herbal-medicines.XLSX"

# interim (tcm)
## latin genitives
### i
pathInterimTcmLatinGenitiveI <-
  "../data/interim/dictionaries/latinGenitiveIDic.tsv"

### is
pathInterimTcmLatinGenitiveIs <-
  "../data/interim/dictionaries/latinGenitiveIsDic.tsv"

## plant parts
pathInterimTcmPlantParts <-
  "../data/interim/dictionaries/tcmPartsManualSubtraction.tsv"

## manual subtraction
pathInterimTcmManualSubtraction <-
  "../data/interim/dictionaries/tcmManualSubtraction.tsv"

## tcm names dic
pathInterimTcmNamesDic <-
  "../data/interim/dictionaries/tcmNamesDic.tsv.zip"

# backbone translation
## manual subtraction
pathInterimTaxaManualSubtraction <-
  "../data/interim/dictionaries/taxaManualSubtraction.tsv"

# interim (backbone)
## problematic names dic
pathInterimProblematicNamesDic <-
  "../data/interim/dictionaries/problematicNamesDic.tsv.zip"

# dbSource
## afrotryp
### original
pathAfrotrypOriginal <-
  "../data/external/dbSource/AFROTRYP/0_initial_files/AFROTRYP.tsv.zip"

### standard
pathAfrotrypStandard <- "../data/interim/db/AFROTRYP_std.tsv.zip"

## alkamid
### original
pathAlkamidOriginal <-
  "../data/external/dbSource/ALKAMID/0_initial_files/ALKAMID_scraped.tsv.zip"

### ref
pathAlkamidRef <-
  "../data/external/dbSource/ALKAMID/0_initial_files/ALKAMID_ref_scraped.tsv.zip"

### standard
pathAlkamidStandard <- "../data/interim/db/ALKAMID_std.tsv.zip"

## biofacquim
### original
pathBiofacquimOriginal <-
  "../data/external/dbSource/BIOFACQUIM/0_initial_files/apps_database_csv_BIOFACQUIM.csv.zip"

### standard
pathBiofacquimStandard <-
  "../data/interim/db/BIOFACQUIM_std.tsv.zip"

## biophytmol
### original
pathBiophytmolOriginal <-
  "../data/external/dbSource/BIOPHYTMOL/0_initial_files/BIOPHYTMOL_scraped.tsv.zip"

### standard
pathBiophytmolStandard <-
  "../data/interim/db/BIOPHYTMOL_std.tsv.zip"

## carotenoiddb
### original
pathCarotenoiddbInchiKey <-
  "../data/external/dbSource/CAROTENOIDDB/0_initial_files/Carotenoids_InChI_InChIKey.tsv.zip"

### original
pathCarotenoiddbOriginal <-
  "../data/external/dbSource/CAROTENOIDDB/0_initial_files/CAROTENOIDDB_scraped.tsv.zip"

### standard
pathCarotenoiddbStandard <-
  "../data/interim/db/CAROTENOIDDB_std.tsv.zip"

## cmaup
### original
#### 1
pathCmaupOriginal_1 <-
  "../data/external/dbSource/CMAUP/0_initial_files/CMAUPv1.0_download_Ingredients_All.txt"

#### 2
pathCmaupOriginal_2 <-
  "../data/external/dbSource/CMAUP/0_initial_files/CMAUPv1.0_download_Plants.txt"

#### 3
pathCmaupOriginal_3 <-
  "../data/external/dbSource/CMAUP/0_initial_files/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt"

### standard
pathCmaupStandard <- "../data/interim/db/CMAUP_std.tsv.zip"

## coconut
### original
pathCoconutOriginal <-
  "../data/external/dbSource/COCONUT/0_initial_files/COCONUT.tsv.zip"

### standard
pathCoconutStandard <- "../data/interim/db/COCONUT_std.tsv.zip"

## cyanometdb
### original
pathCyanometdbOriginal <-
  "../data/external/dbSource/CYANOMETDB/0_initial_files/media-1.csv"

### standard
pathCyanometdbStandard <-
  "../data/interim/db/CYANOMETDB_std.tsv.zip"

## dnp
### original
pathDnpOriginal <-
  "../data/external/dbSource/DNP/0_initial_files/28_2/full_set.csv"

### standard
pathDnpStandard <- "../data/interim/db/DNP_std.tsv.zip"

## drduke
### common names
pathDrdukeCommonNames <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/COMMON_NAMES.csv"

### farmacy
pathDrdukeFarmacy <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/FARMACY_NEW.csv"

### taxa
pathDrdukeTaxa <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/FNFTAX.csv"

### taxa
pathDrdukeRef <-
  "../data/external/dbSource/DRDUKE/0_initial_files/Duke-Source-CSV/REFERENCES.csv"

### standard
pathDrdukeStandard <- "../data/interim/db/DRDUKE_std.tsv.zip"

## etcm
### original
pathEtcmOriginal <-
  list.files(path = "../data/external/dbSource/ETCM/0_initial_files/data/",
             pattern = "*.csv",
             full.names = TRUE)

### standard
pathEtcmStandard <- "../data/interim/db/ETCM_std.tsv.zip"

## foodb
### compounds flavors
pathFoodbCompoundsFlavors <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv"

### compounds
pathFoodbCompounds <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Compound_copy.csv"

### content
pathFoodbContent <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Content.csv"

### flavor
pathFoodbFlavor <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Flavor.csv"

### foods
pathFoodbFoods <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Food_copy.csv"

### ref
pathFoodbRef <-
  "../data/external/dbSource/FOODB/0_initial_files/foodb_2020_04_07_csv/Reference.csv"

### standard
pathFoodbStandard <- "../data/interim/db/FOODB_std.tsv.zip"



### CONTINUE  DBs CLEANING HERE ###



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

## taxa levels dic
pathInterimTaxaLevelsDic <-
  "../data/interim/dictionaries/taxaLevelsDic.tsv"

## black listed string dic
pathInterimBlackDic <-
  "../data/interim/dictionaries/blackDic.tsv"



