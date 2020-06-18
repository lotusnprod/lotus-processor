DATA_PATH ?= ${PWD}/data

INTERIM_PATH = ${DATA_PATH}/interim/db
SOURCE_PATH = ${DATA_PATH}/external/dbSource

AFROTRYP_SOURCE_PATH = ${SOURCE_PATH}/afrotryp
ALKAMID_SOURCE_PATH = ${SOURCE_PATH}/alkamid
BIOFACQUIM_SOURCE_PATH = ${SOURCE_PATH}/biofacquim
BIOPHYTMOL_SOURCE_PATH = ${SOURCE_PATH}/biophytmol
CAROTENOIDDB_SOURCE_PATH = ${SOURCE_PATH}/carotenoiddb
CMAUP_SOURCE_PATH = ${SOURCE_PATH}/cmaup
COCONUT_SOURCE_PATH = ${SOURCE_PATH}/coconut
CYANOMETDB_SOURCE_PATH = ${SOURCE_PATH}/cyanometdb
DNP_SOURCE_PATH = ${SOURCE_PATH}/dnp
DRDUKE_SOURCE_PATH = ${SOURCE_PATH}/drduke
ETCM_SOURCE_PATH = ${SOURCE_PATH}/etcm
FOODB_SOURCE_PATH = ${SOURCE_PATH}/foodb
INFLAMNAT_SOURCE_PATH = ${SOURCE_PATH}/inflamnat
KNAPSACK_SOURCE_PATH = ${SOURCE_PATH}/knapsack
METABOLIGHTS_SOURCE_PATH = ${SOURCE_PATH}/metabolights


.PHONY: help docker-build docker-bash databases afrotryp alkamid alkamid-rescrape biofacquim biophytmol biophytmol-rescrape carotenoiddb carotenoiddb-rescrape cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack knapsack-rescrape metabolights metabolights-rescrape metabolights-reconvert
.PHONY: curating curating-integrating curating-editing curating-editing-bio

help:
	@echo "Builder"
	@echo "-------"
	@echo ""
	@echo "docker-build: build the docker image (with no data)"
	@echo "docker-bash: run a shell into the docker image"
	@echo "databases: build the databases (no scraping)"
	@echo "databases-rescrape: rescrape the databases (when possible)"
	@echo ""
	@echo "curating: Run the 2_curating scripts"

docker-build:
	docker build -t onpdb-environment .

docker-bash:
	docker run -it --rm -v $$PWD:/srv/onpdb onpdb-environment bash

databases: afrotryp alkamid biofacquim biophytmol carotenoiddb cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack metabolights

databases-reconvert: metabolights-reconvert

databases-rescrape: alkamid-rescrape biophytmol-rescrape carotenoiddb-rescrape knapsack-rescrape metabolights-rescrape

afrotryp: ${INTERIM_PATH}/afrotryp.tsv.zip

${INTERIM_PATH}/afrotryp.tsv.zip: ${AFROTRYP_SOURCE_PATH}/afrotryp.tsv.zip
	cd src &&	Rscript 1_gathering/db/afrotryp/standardizing.R

alkamid: ${INTERIM_PATH}/alkamid.tsv.zip

${INTERIM_PATH}/alkamid.tsv.zip: ${ALKAMID_SOURCE_PATH}/alkamidRefScraped.tsv.zip ${ALKAMID_SOURCE_PATH}/alkamidScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/alkamid/standardizing.R

alkamid-rescrape:
	cd src && Rscript 1_gathering/db/alkamid/scraping.R

biofacquim: ${INTERIM_PATH}/biofacquim.tsv.zip

${INTERIM_PATH}/biofacquim.tsv.zip: ${BIOFACQUIM_SOURCE_PATH}/apps_database_csv_BIOFACQUIM.csv
	cd src &&	Rscript 1_gathering/db/biofacquim/standardizing.R

biophytmol: ${INTERIM_PATH}/biophytmol.tsv.zip

${INTERIM_PATH}/biophytmol.tsv.zip: ${BIOPHYTMOL_SOURCE_PATH}/biophytmolScraped.tsv.zip 
	cd src &&	Rscript 1_gathering/db/biophytmol/standardizing.R

biophytmol-rescrape:
	cd src && Rscript 1_gathering/db/biophytmol/scraping.R

carotenoiddb: ${INTERIM_PATH}/carotenoiddb.tsv.zip

${INTERIM_PATH}/carotenoiddb.tsv.zip: ${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip 
	cd src &&	Rscript 1_gathering/db/carotenoiddb/standardizing.R

carotenoiddb-rescrape: ${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip 

${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip : ${CAROTENOIDDB_SOURCE_PATH}/Carotenoids_InChI_InChIKey.tsv
	cd src && Rscript 1_gathering/db/carotenoiddb/scraping.R

cmaup: ${INTERIM_PATH}/cmaup.tsv.zip

${INTERIM_PATH}/cmaup.tsv.zip: ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Ingredients_All.txt ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Plants.txt ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt
	cd src &&	Rscript 1_gathering/db/cmaup/standardizing.R

# maybe not the right way to do it
coconut: ${INTERIM_PATH}/coconut.tsv.zip

${INTERIM_PATH}/coconut.tsv.zip: ${COCONUT_SOURCE_PATH}/COCONUT.sdf.zip ${COCONUT_SOURCE_PATH}/coconutConverted.tsv.zip
	cd src &&	python 1_gathering/db/coconut/converting.py  && Rscript 1_gathering/db/coconut/standardizing.R

cyanometdb: ${INTERIM_PATH}/cyanometdb.tsv.zip

${INTERIM_PATH}/cyanometdb.tsv.zip: ${CYANOMETDB_SOURCE_PATH}/media-1.csv
	cd src &&	Rscript 1_gathering/db/cyanometdb/standardizing.R

dnp: ${INTERIM_PATH}/dnp.tsv.zip

${INTERIM_PATH}/dnp.tsv.zip: ${DNP_SOURCE_PATH}/28_2/full_set.csv
	cd src &&	Rscript 1_gathering/db/dnp/standardizing.R

drduke: ${INTERIM_PATH}/drduke.tsv.zip

${INTERIM_PATH}/drduke.tsv.zip: ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/COMMON_NAMES.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/FARMACY_NEW.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/FNFTAX.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/REFERENCES.csv
	cd src &&	Rscript 1_gathering/db/drduke/standardizing.R

etcm: ${INTERIM_PATH}/etcm.tsv.zip

${INTERIM_PATH}/etcm.tsv.zip: $(wildcard ${ETCM_SOURCE_PATH}data/*.csv)
	cd src &&	Rscript 1_gathering/db/etcm/standardizing.R

foodb: ${INTERIM_PATH}/foodb.tsv.zip

${INTERIM_PATH}/foodb.tsv.zip: ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Compound_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Content.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Flavor.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Food_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Reference.csv
	cd src &&	Rscript 1_gathering/db/foodb/standardizing.R

inflamnat: ${INTERIM_PATH}/inflamnat.tsv.zip

${INTERIM_PATH}/inflamnat.tsv.zip: ${INFLAMNAT_SOURCE_PATH}/ci8b00560_si_001.xlsx
	cd src &&	Rscript 1_gathering/db/inflamnat/standardizing.R

knapsack: ${INTERIM_PATH}/knapsack.tsv.zip

${INTERIM_PATH}/knapsack.tsv.zip: ${KNAPSACK_SOURCE_PATH}/knapsackScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/knapsack/standardizing.R

knapsack-rescrape:
	cd src && Rscript 1_gathering/db/knapsack/scraping.R


# PROCESS DONE FOR METABOLIGHTS IS HORROR... my apologies I can't write it another way
metabolights: ${INTERIM_PATH}/metabolights.tsv.zip

${INTERIM_PATH}/metabolights.tsv.zip: ${METABOLIGHTS_SOURCE_PATH}/metabolightsPrecleaned.tsv.zip ${METABOLIGHTS_SOURCE_PATH}/metabolightsStudiesScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/metabolights/standardizing.R

metabolights-reconvert: ${INTERIM_PATH}/metabolightsPrecleaned.tsv.zip ${METABOLIGHTS_SOURCE_PATH}/metabolightsStudiesScraped.tsv.zip

${INTERIM_PATH}/metabolightsPrecleaned.tsv.zip: $(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json) ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_complete.xml
	cd src &&	Rscript 1_gathering/db/metabolights/prestandardizing.R &&	Rscript 1_gathering/db/metabolights/standardizingStudies.R

metabolights-rescrape: $(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json)

$(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json) : ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_studies.xml
	cd src && Rscript 1_gathering/db/metabolights/scraping.R


curating: curating-integrating curating-editing

curating-integrating:
	cd src && Rscript 2_curating/1_integrating/integratingOriginalDatabase.R

curating-editing:  curating-editing-bio

curating-editing-bio:
	cd src && Rscript 2_curating/2_editing/bio/editing.R
