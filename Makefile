DATA_PATH ?= ${PWD}/data
SRC_PATH ?= ${PWD}/src

INTERIM_PATH = ${DATA_PATH}/interim/db
SOURCE_PATH = ${DATA_PATH}/external/dbSource
# SOURCE_PATH is somewhat ambiguous .... I was expecting the src directory here ...
## be carfeul with comments they can ruin the makefile as stated here: https://www.gnu.org/software/make/manual/html_node/Recipe-Syntax.html. At least it made mine do not work anymore

# Below is and endless PATH nightmare. Why dont we have a unique ITERIM folder to dump all files, uniquely named ?
## This is because the files should be slowly removed and renamed in dictionaries instead. so working with interim/original table at the begining, bunch of dictionaries (original to translated, translated to curated etc) and processed table in the end, in my view at least
INTERIM_TABLE_PATH = ${DATA_PATH}/interim/tables
INTERIM_TABLE_ORIGINAL_PATH = ${INTERIM_TABLE_PATH}/0_original
INTERIM_TABLE_TRANSLATED_PATH = ${INTERIM_TABLE_PATH}/1_translated
INTERIM_TABLE_CLEANED_PATH = ${INTERIM_TABLE_PATH}/2_cleaned
INTERIM_TABLE_CURATED_PATH = ${INTERIM_TABLE_PATH}/3_curated

SRC_GATHERING_PATH = ${SRC_PATH}/1_gathering/db
SRC_CURATING_PATH = ${SRC_PATH}/2_curating

SRC_CURATING_EDITING_PATH = ${SRC_CURATING_PATH}/2_editing
SRC_CURATING_EDITING_CHEMO_PATH =  ${SRC_CURATING_EDITING_PATH}/chemo
SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_PATH =  ${SRC_CURATING_EDITING_CHEMO_PATH}/subscripts
SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH =  ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_PATH}/1_translating

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

.PHONY: help docker-build docker-bash databases afrotryp alkamid alkamid-rescrape biofacquim biofacquim-reconvert biophytmol biophytmol-rescrape carotenoiddb carotenoiddb-rescrape cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack knapsack-rescrape metabolights metabolights-rescrape metabolights-reconvert mibig mitishamba mitishamba-rescrape nanpdb nanpdb-rescrape npass npatlas npcare npedia npedia-rescrape pamdb phenolexplorer phytohub phytohub-rescrape plantcyc plantcyc-reintegrate procardb procardb-rescrape respect sancdb sancdb-rescrape streptomedb streptomedb-reconvert swmd swmd-rescrape symmap tmdb tmdb-rescrape tmmc tppt triforc triforc-reintegrate unpd unpd-reintegrate
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

databases-reconvert: biofacquim-reconvert metabolights-reconvert streptomedb-reconvert

databases-reintegrate: plantcyc-reintegrate triforc-reintegrate unpd-reintegrate

databases-rescrape: alkamid-rescrape biophytmol-rescrape carotenoiddb-rescrape knapsack-rescrape metabolights-rescrape mitishamba-rescrape nanpdb-rescrape npedia-rescrape phytohub-rescrape procardb-rescrape sancdb-rescrape swmd-rescrape tmdb-rescrape

afrotryp: ${INTERIM_PATH}/afrotryp.tsv.zip

${INTERIM_PATH}/afrotryp.tsv.zip: ${AFROTRYP_SOURCE_PATH}/afrotryp.tsv.zip ${SRC_GATHERING_PATH}/afrotryp/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/afrotryp/standardizing.R

alkamid: ${INTERIM_PATH}/alkamid.tsv.zip

${INTERIM_PATH}/alkamid.tsv.zip: ${ALKAMID_SOURCE_PATH}/alkamidRefScraped.tsv.zip ${ALKAMID_SOURCE_PATH}/alkamidScraped.tsv.zip ${SRC_GATHERING_PATH}/alkamid/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/alkamid/standardizing.R

alkamid-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/alkamid/scraping.R

biofacquim: ${INTERIM_PATH}/biofacquim.tsv.zip

${INTERIM_PATH}/biofacquim.tsv.zip: ${BIOFACQUIM_SOURCE_PATH}/biofacquim.tsv.zip ${SRC_GATHERING_PATH}/biofacquim/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/biofacquim/standardizing.R

biofacquim-reconvert: ${BIOFACQUIM_SOURCE_PATH}/biofacquim.tsv.zip

${BIOFACQUIM_SOURCE_PATH}/biofacquim.tsv.zip: ${BIOFACQUIM_SOURCE_PATH}/BIOFACQUIM_V2.sdf ${SRC_GATHERING_PATH}/biofacquim/converting.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/biofacquim/converting.R

biophytmol: ${INTERIM_PATH}/biophytmol.tsv.zip

${INTERIM_PATH}/biophytmol.tsv.zip: ${BIOPHYTMOL_SOURCE_PATH}/biophytmolScraped.tsv.zip ${SRC_GATHERING_PATH}/biophytmol/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/biophytmol/standardizing.R

biophytmol-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/biophytmol/scraping.R

carotenoiddb: ${INTERIM_PATH}/carotenoiddb.tsv.zip

${INTERIM_PATH}/carotenoiddb.tsv.zip: ${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip ${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R

carotenoiddb-rescrape: ${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip

${CAROTENOIDDB_SOURCE_PATH}/carotenoiddbScraped.tsv.zip : ${CAROTENOIDDB_SOURCE_PATH}/Carotenoids_InChI_InChIKey.tsv ${SRC_GATHERING_PATH}/carotenoiddb/scraping.R
	cd src && Rscript ${SRC_GATHERING_PATH}/carotenoiddb/scraping.R

cmaup: ${INTERIM_PATH}/cmaup.tsv.zip

${INTERIM_PATH}/cmaup.tsv.zip: ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Ingredients_All.txt ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Plants.txt ${CMAUP_SOURCE_PATH}/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt ${SRC_GATHERING_PATH}/cmaup/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/cmaup/standardizing.R

# maybe not the right way to do it
coconut: ${INTERIM_PATH}/coconut.tsv.zip

${INTERIM_PATH}/coconut.tsv.zip: ${COCONUT_SOURCE_PATH}/COCONUT.sdf.zip ${COCONUT_SOURCE_PATH}/coconutConverted.tsv.zip ${SRC_GATHERING_PATH}/coconut/converting.py ${SRC_GATHERING_PATH}/coconut/standardizing.R
	cd src &&	python ${SRC_GATHERING_PATH}/coconut/converting.py  && Rscript ${SRC_GATHERING_PATH}/coconut/standardizing.R

cyanometdb: ${INTERIM_PATH}/cyanometdb.tsv.zip

${INTERIM_PATH}/cyanometdb.tsv.zip: ${CYANOMETDB_SOURCE_PATH}/media-1.csv ${SRC_GATHERING_PATH}/cyanometdb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/cyanometdb/standardizing.R

dnp: ${INTERIM_PATH}/dnp.tsv.zip

${INTERIM_PATH}/dnp.tsv.zip: ${DNP_SOURCE_PATH}/28_2/full_set.csv ${SRC_GATHERING_PATH}/dnp/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/dnp/standardizing.R

drduke: ${INTERIM_PATH}/drduke.tsv.zip

${INTERIM_PATH}/drduke.tsv.zip: ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/COMMON_NAMES.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/FARMACY_NEW.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/FNFTAX.csv ${DRDUKE_SOURCE_PATH}/Duke-Source-CSV/REFERENCES.csv ${SRC_GATHERING_PATH}/drduke/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/drduke/standardizing.R

etcm: ${INTERIM_PATH}/etcm.tsv.zip

${INTERIM_PATH}/etcm.tsv.zip: $(wildcard ${ETCM_SOURCE_PATH}/data/*.csv) ${SRC_GATHERING_PATH}/etcm/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/etcm/standardizing.R

foodb: ${INTERIM_PATH}/foodb.tsv.zip

${INTERIM_PATH}/foodb.tsv.zip: ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Compound_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Content.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Flavor.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Food_copy.csv ${FOODB_SOURCE_PATH}/foodb_2020_04_07_csv/Reference.csv ${SRC_GATHERING_PATH}/foodb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/foodb/standardizing.R

inflamnat: ${INTERIM_PATH}/inflamnat.tsv.zip

${INTERIM_PATH}/inflamnat.tsv.zip: ${INFLAMNAT_SOURCE_PATH}/ci8b00560_si_001.xlsx ${SRC_GATHERING_PATH}/inflamnat/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/inflamnat/standardizing.R

knapsack: ${INTERIM_PATH}/knapsack.tsv.zip

${INTERIM_PATH}/knapsack.tsv.zip: ${KNAPSACK_SOURCE_PATH}/knapsackScraped.tsv.zip ${SRC_GATHERING_PATH}/knapsack/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/knapsack/standardizing.R

knapsack-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/knapsack/scraping.R

# PROCESS DONE FOR METABOLIGHTS IS HORROR... my apologies I can't write it another way
metabolights: ${INTERIM_PATH}/metabolights.tsv.zip

${INTERIM_PATH}/metabolights.tsv.zip: ${METABOLIGHTS_SOURCE_PATH}/metabolightsPrecleaned.tsv.zip ${METABOLIGHTS_SOURCE_PATH}/metabolightsStudiesScraped.tsv.zip ${SRC_GATHERING_PATH}/metabolights/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/metabolights/standardizing.R

metabolights-reconvert: ${METABOLIGHTS_SOURCE_PATH}/metabolightsPrecleaned.tsv.zip ${METABOLIGHTS_SOURCE_PATH}/metabolightsStudiesScraped.tsv.zip

${METABOLIGHTS_SOURCE_PATH}/metabolightsPrecleaned.tsv.zip: $(wildcard ${METABOLIGHTS_SOURCE_PATH}/studiesScraped/*.json) ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_complete.xml ${SRC_GATHERING_PATH}/metabolights/prestandardizing.R ${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/metabolights/prestandardizing.R &&	Rscript ${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R

metabolights-rescrape: $(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json)

$(wildcard ${METABOLIGHTS_SOURCE_PATH}studiesScraped/*.json) : ${METABOLIGHTS_SOURCE_PATH}/eb-eye_metabolights_studies.xml ${SRC_GATHERING_PATH}/metabolights/scraping.R
	cd src && Rscript ${SRC_GATHERING_PATH}/metabolights/scraping.R

mibig: ${INTERIM_PATH}/mibig.tsv.zip

${INTERIM_PATH}/mibig.tsv.zip: $(wildcard ${MIBIG_SOURCE_PATH}/mibig_json_2.0/*.json) ${SRC_GATHERING_PATH}/mibig/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/mibig/standardizing.R

mitishamba: ${INTERIM_PATH}/mitishamba.tsv.zip

${INTERIM_PATH}/mitishamba.tsv.zip: ${MITISHAMBA_SOURCE_PATH}/mitishambaScraped.tsv.zip ${SRC_GATHERING_PATH}/mitishamba/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/mitishamba/standardizing.R

mitishamba-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/mitishamba/scraping.R

nanpdb: ${INTERIM_PATH}/nanpdb.tsv.zip

${INTERIM_PATH}/nanpdb.tsv.zip: ${NANPDB_SOURCE_PATH}/nanpdbScraped.tsv.zip ${SRC_GATHERING_PATH}/nanpdb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/nanpdb/standardizing.R

nanpdb-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/nanpdb/scraping.R

npass: ${INTERIM_PATH}/npass.tsv.zip

${INTERIM_PATH}/npass.tsv.zip: ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_generalInfo.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_properties.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_speciesInfo.txt ${NPASS_SOURCE_PATH}/NPASSv1.0_download_naturalProducts_species_pair.txt ${SRC_GATHERING_PATH}/npass/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/npass/standardizing.R

npatlas: ${INTERIM_PATH}/npatlas.tsv.zip

${INTERIM_PATH}/npatlas.tsv.zip: ${NPATLAS_SOURCE_PATH}/np_atlas_2019_12.tsv ${SRC_GATHERING_PATH}/npatlas/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/npatlas/standardizing.R

npcare: ${INTERIM_PATH}/npcare.tsv.zip

${INTERIM_PATH}/npcare.tsv.zip: ${NPCARE_SOURCE_PATH}/npcare.zip ${SRC_GATHERING_PATH}/npcare/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/npcare/standardizing.R

npedia: ${INTERIM_PATH}/npedia.tsv.zip

${INTERIM_PATH}/npedia.tsv.zip: ${NPEDIA_SOURCE_PATH}/npediaScraped.tsv.zip ${SRC_GATHERING_PATH}/npedia/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/npedia/standardizing.R

npedia-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/npedia/scraping.R

pamdb: ${INTERIM_PATH}/pamdb.tsv.zip

${INTERIM_PATH}/pamdb.tsv.zip: ${PAMDB_SOURCE_PATH}/PaMet.xlsx ${SRC_GATHERING_PATH}/pamdb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/pamdb/standardizing.R

phenolexplorer: ${INTERIM_PATH}/phenolexplorer.tsv.zip

${INTERIM_PATH}/phenolexplorer.tsv.zip: ${PHENOLEXPLORER_SOURCE_PATH}/compounds-classification.csv ${PHENOLEXPLORER_SOURCE_PATH}/compounds-structures.csv ${PHENOLEXPLORER_SOURCE_PATH}/compounds.csv ${PHENOLEXPLORER_SOURCE_PATH}/foods-classification.csv ${PHENOLEXPLORER_SOURCE_PATH}/foods.csv ${PHENOLEXPLORER_SOURCE_PATH}/metabolites-structures.csv ${PHENOLEXPLORER_SOURCE_PATH}/metabolites.csv ${PHENOLEXPLORER_SOURCE_PATH}/publications.csv ${PHENOLEXPLORER_SOURCE_PATH}/composition-data.xlsx ${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R

phytohub: ${INTERIM_PATH}/phytohub.tsv.zip

${INTERIM_PATH}/phytohub.tsv.zip: ${PHYTOHUB_SOURCE_PATH}/phytohubScraped.tsv.zip ${SRC_GATHERING_PATH}/phytohub/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/phytohub/standardizing.R

phytohub-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/phytohub/scraping.R

plantcyc: ${INTERIM_PATH}/plantcyc.tsv.zip

${INTERIM_PATH}/plantcyc.tsv.zip: $(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip) ${SRC_GATHERING_PATH}/plantcyc/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/plantcyc/standardizing.R

# do not know how to proceed here... might be wrong
plantcyc-reintegrate: $(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip)

$(wildcard ${PLANTCYC_SOURCE_PATH}/*.tsv.zip) : $(wildcard ${PLANTCYC_SOURCE_PATH}/0_data/*/*/data/compounds.dat) ${SRC_GATHERING_PATH}/plantcyc/integrating.R
	cd src && Rscript ${SRC_GATHERING_PATH}/plantcyc/integrating.R

procardb: ${INTERIM_PATH}/procardb.tsv.zip

${INTERIM_PATH}/procardb.tsv.zip: ${PROCARDB_SOURCE_PATH}/procardbScraped.tsv.zip ${SRC_GATHERING_PATH}/procardb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/procardb/standardizing.R

prodardb-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/procardb/scraping.R

respect: ${INTERIM_PATH}/respect.tsv.zip

${INTERIM_PATH}/respect.tsv.zip: $(wildcard ${RESPECT_SOURCE_PATH}/respect/*.txt) ${SRC_GATHERING_PATH}/respect/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/respect/standardizing.R

sancdb: ${INTERIM_PATH}/sancdb.tsv.zip

${INTERIM_PATH}/sancdb.tsv.zip: ${SANCDB_SOURCE_PATH}/sancdbScraped.tsv.zip ${SRC_GATHERING_PATH}/sancdb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/sancdb/standardizing.R

sancdb-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/sancdb/scraping.R

streptomedb: ${INTERIM_PATH}/streptomedb.tsv.zip

${INTERIM_PATH}/streptomedb.tsv.zip: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip ${SRC_GATHERING_PATH}/streptomedb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/streptomedb/standardizing.R

streptomedb-reconvert: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip

${STREPTOMEDB_SOURCE_PATH}/streptomedb.tsv.zip: ${STREPTOMEDB_SOURCE_PATH}/streptomedb.sdf ${SRC_GATHERING_PATH}/streptomedb/converting.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/streptomedb/converting.R

swmd: ${INTERIM_PATH}/swmd.tsv.zip

${INTERIM_PATH}/swmd.tsv.zip: ${SWMD_SOURCE_PATH}/swmdScraped.tsv.zip ${SRC_GATHERING_PATH}/swmd/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/swmd/standardizing.R

swmd-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/swmd/scraping.R

symmap: ${INTERIM_PATH}/symmap.tsv.zip

${INTERIM_PATH}/symmap.tsv.zip: $(wildcard ${SYMMAP_SOURCE_PATH}/data/*.csv) ${SRC_GATHERING_PATH}/symmap/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/symmap/standardizing.R

tmdb: ${INTERIM_PATH}/tmdb.tsv.zip

${INTERIM_PATH}/tmdb.tsv.zip: ${TMDB_SOURCE_PATH}/tmdbScraped.tsv.zip ${SRC_GATHERING_PATH}/tmdb/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/tmdb/standardizing.R

tmdb-rescrape:
	cd src && Rscript ${SRC_GATHERING_PATH}/tmdb/scraping.R

tmmc: ${INTERIM_PATH}/tmmc.tsv.zip

${INTERIM_PATH}/tmmc.tsv.zip: ${TMMC_SOURCE_PATH}/compound.xlsx ${SRC_GATHERING_PATH}/tmmc/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/tmmc/standardizing.R

tppt: ${INTERIM_PATH}/tppt.tsv.zip

${INTERIM_PATH}/tppt.tsv.zip: ${TPPT_SOURCE_PATH}/TPPT_database.xlsx ${SRC_GATHERING_PATH}/tppt/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/tppt/standardizing.R

triforc: ${INTERIM_PATH}/triforc.tsv.zip

${INTERIM_PATH}/triforc.tsv.zip: ${TRIFORC_SOURCE_PATH}/triforcBis.tsv ${SRC_GATHERING_PATH}/triforc/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/triforc/standardizing.R

triforc-reintegrate: ${TRIFORC_SOURCE_PATH}/triforcBis.tsv

${TRIFORC_SOURCE_PATH}/triforcBis.tsv : ${TRIFORC_SOURCE_PATH}/triforcOriginal.tsv ${TRIFORC_SOURCE_PATH}/triforcToGet.tsv ${SRC_GATHERING_PATH}/triforc/prestandardizing.R
	cd src && Rscript ${SRC_GATHERING_PATH}/triforc/prestandardizing.R

unpd: ${INTERIM_PATH}/unpd.tsv.zip

${INTERIM_PATH}/unpd.tsv.zip: ${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip ${SRC_GATHERING_PATH}/unpd/standardizing.R
	cd src &&	Rscript ${SRC_GATHERING_PATH}/unpd/standardizing.R

unpd-reintegrate: ${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip

${UNPD_SOURCE_PATH}/unpdIntegrated.tsv.zip: ${UNPD_SOURCE_PATH}/unpd_final.csv.zip ${UNPD_SOURCE_PATH}/UNPD_DB.csv.zip ${SRC_GATHERING_PATH}/unpd/integrating.R
	cd src && Rscript ${SRC_GATHERING_PATH}/unpd/integrating.R

curating: curating-integrating curating-editing

curating-integrating:
	cd src && Rscript 2_curating/1_integrating/integratingOriginalDatabase.R

curating-editing: curating-editing-bio curating-editing-reference

curating-editing-bio:
	cd src && Rscript 2_curating/2_editing/bio/editing.R

curating-editing-reference:
	cd src && Rscript 2_curating/2_editing/reference/editing.R

curating-editing-chemo: curating-editing-chemo-name curating-editing-chemo-smiles

curating-editing-chemo-name: ${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip

${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureNominal.tsv.zip: ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R
${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/names.R: ${INTERIM_TABLE_ORIGINAL_PATH}/originalStructureNominal.tsv.zip
	cd src && Rscript 2_curating/2_editing/chemo/subscripts/1_translating/names.R

curating-editing-chemo-smiles: 
	cd src && python 2_curating/2_editing/chemo/subscripts/1_translating/smiles.py
