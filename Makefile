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
MIBIG_SOURCE_PATH = ${SOURCE_PATH}/mibig
MITISHAMBA_SOURCE_PATH = ${SOURCE_PATH}/mitishamba
NANPDB_SOURCE_PATH = ${SOURCE_PATH}/nanpdb
NPASS_SOURCE_PATH = ${SOURCE_PATH}/npass
NPATLAS_SOURCE_PATH = ${SOURCE_PATH}/npatlas
NPCARE_SOURCE_PATH = ${SOURCE_PATH}/npcare
NPEDIA_SOURCE_PATH = ${SOURCE_PATH}/npedia
PAMDB_SOURCE_PATH = ${SOURCE_PATH}/pamdb
PHENOLEXPLORER_SOURCE_PATH = ${SOURCE_PATH}/phenolexplorer
PHYTOHUB_SOURCE_PATH = ${SOURCE_PATH}/phytohub
PLANTCYC_SOURCE_PATH = ${SOURCE_PATH}/plantcyc
PROCARDB_SOURCE_PATH = ${SOURCE_PATH}/procardb
RESPECT_SOURCE_PATH = ${SOURCE_PATH}/respect
SANCDB_SOURCE_PATH = ${SOURCE_PATH}/sancdb
STREPTOMEDB_SOURCE_PATH = ${SOURCE_PATH}/streptomedb
SWMD_SOURCE_PATH = ${SOURCE_PATH}/swmd
SYMMAP_SOURCE_PATH = ${SOURCE_PATH}/symmap
TMDB_SOURCE_PATH = ${SOURCE_PATH}/tmdb
TMMC_SOURCE_PATH = ${SOURCE_PATH}/tmmc
TPPT_SOURCE_PATH = ${SOURCE_PATH}/tppt
TRIFORC_SOURCE_PATH = ${SOURCE_PATH}/triforc
UNPD_SOURCE_PATH = ${SOURCE_PATH}/unpd


.PHONY: help docker-build docker-bash databases afrotryp alkamid alkamid-rescrape biofacquim biophytmol biophytmol-rescrape carotenoiddb carotenoiddb-rescrape cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack knapsack-rescrape metabolights metabolights-rescrape metabolights-reconvert mibig mitishamba mitishamba-rescrape nanpdb nanpdb-rescrape npass npatlas npcare npedia npedia-rescrape pamdb phenolexplorer phytohub phytohub-rescrape plantcyc plantcyc-reintegrate procardb procardb-rescrape respect sancdb sancdb-rescrape streptomedb streptomedb-reconvert swmd swmd-rescrape symmap tmdb tmdb-rescrape tmmc tppt triforc triforc-reintegrate unpd unpd-reintegrate
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

databases: afrotryp alkamid biofacquim biophytmol carotenoiddb cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack metabolights mibig mitishamba nanpdb npass npatlas npcare npedia pamdb phenolexplorer phytohub plantcyc procardb respect sancdb streptomedb swmd symmap tmdb tmmc tppt triforc unpd

databases-reconvert: metabolights-reconvert streptomedb-reconvert

databases-reintegrate: plantcyc-reintegrate triforc-reintegrate unpd-reintegrate

databases-rescrape: alkamid-rescrape biophytmol-rescrape carotenoiddb-rescrape knapsack-rescrape metabolights-rescrape mitishamba-rescrape nanpdb-rescrape npedia-rescrape phytohub-rescrape procardb-rescrape sancdb-rescrape swmd-rescrape tmdb-rescrape

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

${INTERIM_PATH}/etcm.tsv.zip: $(wildcard ${ETCM_SOURCE_PATH}/data/*.csv)
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

${INTERIM_PATH}/metabolightsPrecleaned.tsv.zip: $(wildcard ${METABOLIGHTS_SOURCE_PATH}/studiesScraped/*.json) ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_complete.xml
	cd src &&	Rscript 1_gathering/db/metabolights/prestandardizing.R &&	Rscript 1_gathering/db/metabolights/standardizingStudies.R

metabolights-rescrape: $(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json)

$(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json) : ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_studies.xml
	cd src && Rscript 1_gathering/db/metabolights/scraping.R

mibig: ${INTERIM_PATH}/mibig.tsv.zip

${INTERIM_PATH}/mibig.tsv.zip: $(wildcard ${MIBIG_SOURCE_PATH}/mibig_json_2.0/*.json)
	cd src &&	Rscript 1_gathering/db/mibig/standardizing.R

mitishamba: ${INTERIM_PATH}/mitishamba.tsv.zip

${INTERIM_PATH}/mitishamba.tsv.zip: ${MITISHAMBA_SOURCE_PATH}/mitishambaScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/mitishamba/standardizing.R

mitishamba-rescrape:
	cd src && Rscript 1_gathering/db/mitishamba/scraping.R

nanpdb: ${INTERIM_PATH}/nanpdb.tsv.zip

${INTERIM_PATH}/nanpdb.tsv.zip: ${NANPDB_SOURCE_PATH}/nanpdbScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/nanpdb/standardizing.R

nanpdb-rescrape:
	cd src && Rscript 1_gathering/db/nanpdb/scraping.R

npass: ${INTERIM_PATH}/npass.tsv.zip

${INTERIM_PATH}/npass.tsv.zip: ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_generalInfo.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_properties.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_speciesInfo.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_species_pair.txt
	cd src &&	Rscript 1_gathering/db/npass/standardizing.R

npatlas: ${INTERIM_PATH}/npatlas.tsv.zip

${INTERIM_PATH}/npatlas.tsv.zip: ${NPATLAS_SOURCE_PATH}/np_atlas_2019_12.tsv
	cd src &&	Rscript 1_gathering/db/npatlas/standardizing.R

npcare: ${INTERIM_PATH}/npcare.tsv.zip

${INTERIM_PATH}/npcare.tsv.zip: ${NPCARE_SOURCE_PATH}/npcare.zip
	cd src &&	Rscript 1_gathering/db/npcare/standardizing.R

npedia: ${INTERIM_PATH}/npedia.tsv.zip

${INTERIM_PATH}/npedia.tsv.zip: ${NPEDIA_SOURCE_PATH}/npediaScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/npedia/standardizing.R

npedia-rescrape:
	cd src && Rscript 1_gathering/db/npedia/scraping.R

pamdb: ${INTERIM_PATH}/pamdb.tsv.zip

${INTERIM_PATH}/pamdb.tsv.zip: ${PAMDB_SOURCE_PATH}/PaMet.xlsx
	cd src &&	Rscript 1_gathering/db/pamdb/standardizing.R

phenolexplorer: ${INTERIM_PATH}/phenolexplorer.tsv.zip

${INTERIM_PATH}/phenolexplorer.tsv.zip: ${PHENOLEXPLORER_SOURCE_PATH}/compounds-classification.csv ${PHENOLEXPLORER_SOURCE_PATH}/compounds-structures.csv ${PHENOLEXPLORER_SOURCE_PATH}/compounds.csv ${PHENOLEXPLORER_SOURCE_PATH}/foods-classification.csv ${PHENOLEXPLORER_SOURCE_PATH}/foods.csv ${PHENOLEXPLORER_SOURCE_PATH}/metabolites-structures.csv ${PHENOLEXPLORER_SOURCE_PATH}/metabolites.csv ${PHENOLEXPLORER_SOURCE_PATH}/publications.csv ${PHENOLEXPLORER_SOURCE_PATH}/composition-data.xlsx
	cd src &&	Rscript 1_gathering/db/phenolexplorer/standardizing.R

phytohub: ${INTERIM_PATH}/phytohub.tsv.zip

${INTERIM_PATH}/phytohub.tsv.zip: ${PHYTOHUB_SOURCE_PATH}/phytohubScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/phytohub/standardizing.R

phytohub-rescrape:
	cd src && Rscript 1_gathering/db/phytohub/scraping.R

plantcyc: ${INTERIM_PATH}/plantcyc.tsv.zip

${INTERIM_PATH}/plantcyc.tsv.zip: $(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip)
	cd src &&	Rscript 1_gathering/db/plantcyc/standardizing.R

# do not know how to proceed here... might be wrong
plantcyc-reintegrate: $(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip)

$(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip) : $(wildcard ${PLANTCYC_SOURCE_PATH}/0_data/*/*/data/compounds.dat)
	cd src && Rscript 1_gathering/db/plantcyc/integrating.R

procardb: ${INTERIM_PATH}/procardb.tsv.zip

${INTERIM_PATH}/procardb.tsv.zip: ${PROCARDB_SOURCE_PATH}/procardbScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/procardb/standardizing.R

prodardb-rescrape:
	cd src && Rscript 1_gathering/db/procardb/scraping.R

respect: ${INTERIM_PATH}/respect.tsv.zip

${INTERIM_PATH}/respect.tsv.zip: $(wildcard ${RESPECT_SOURCE_PATH}/respect/*.txt)
	cd src &&	Rscript 1_gathering/db/respect/standardizing.R

sancdb: ${INTERIM_PATH}/sancdb.tsv.zip

${INTERIM_PATH}/sancdb.tsv.zip: ${SANCDB_SOURCE_PATH}/sancdbScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/sancdb/standardizing.R

sancdb-rescrape:
	cd src && Rscript 1_gathering/db/sancdb/scraping.R

streptomedb: ${INTERIM_PATH}/streptomedb.tsv.zip

${INTERIM_PATH}/streptomedb.tsv.zip: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip
	cd src &&	Rscript 1_gathering/db/streptomedb/standardizing.R

streptomedb-reconvert: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip

${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.sdf
	cd src &&	Rscript 1_gathering/db/streptomedb/converting.R

swmd: ${INTERIM_PATH}/swmd.tsv.zip

${INTERIM_PATH}/swmd.tsv.zip: ${SWMD_SOURCE_PATH}/swmdScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/swmd/standardizing.R

swmd-rescrape:
	cd src && Rscript 1_gathering/db/swmd/scraping.R

symmap: ${INTERIM_PATH}/symmap.tsv.zip

${INTERIM_PATH}/symmap.tsv.zip: $(wildcard ${SYMMAP_SOURCE_PATH}/data/*.csv)
	cd src &&	Rscript 1_gathering/db/symmap/standardizing.R

tmdb: ${INTERIM_PATH}/tmdb.tsv.zip

${INTERIM_PATH}/tmdb.tsv.zip: ${TMDB_SOURCE_PATH}/tmdbScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/tmdb/standardizing.R

tmdb-rescrape:
	cd src && Rscript 1_gathering/db/tmdb/scraping.R

tmmc: ${INTERIM_PATH}/tmmc.tsv.zip

${INTERIM_PATH}/tmmc.tsv.zip: ${TMMC_SOURCE_PATH}/compound.xlsx
	cd src &&	Rscript 1_gathering/db/tmmc/standardizing.R

tppt: ${INTERIM_PATH}/tppt.tsv.zip

${INTERIM_PATH}/tppt.tsv.zip: ${TPPT_SOURCE_PATH}/TPPT_database.xlsx
	cd src &&	Rscript 1_gathering/db/tppt/standardizing.R

triforc: ${INTERIM_PATH}/triforc.tsv.zip

${INTERIM_PATH}/triforc.tsv.zip: ${TRIFORC_SOURCE_PATH}/triforcBis.tsv
	cd src &&	Rscript 1_gathering/db/triforc/standardizing.R

triforc-reintegrate: ${TRIFORC_SOURCE_PATH}/triforcBis.tsv

${TRIFORC_SOURCE_PATH}/triforcBis.tsv : ${TRIFORC_SOURCE_PATH}/triforcOriginal.tsv ${TRIFORC_SOURCE_PATH}/triforcToGet.tsv
	cd src && Rscript 1_gathering/db/triforc/prestandardizing.R

unpd: ${INTERIM_PATH}/unpd.tsv.zip

${INTERIM_PATH}/unpd.tsv.zip: ${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip
	cd src &&	Rscript 1_gathering/db/unpd/standardizing.R

triforc-reintegrate: ${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip
UNPD_SOURCE_PATH
${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip: ${UNPD_SOURCE_PATH}/unpd_final.csv.zip ${UNPD_SOURCE_PATH}/UNPD_DB.csv.zip
	cd src && Rscript 1_gathering/db/unpd/integrating.R	

curating: curating-integrating curating-editing

curating-integrating:
	cd src && Rscript 2_curating/1_integrating/integratingOriginalDatabase.R

curating-editing:  curating-editing-bio

curating-editing-bio:
	cd src && Rscript 2_curating/2_editing/bio/editing.R
