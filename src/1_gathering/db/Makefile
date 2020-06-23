include ../../../paths.mk

DATABASES = afrotryp alkamid biofacquim biophytmol carotenoiddb cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack mibig mitishamba nanpdb npass npatlas npcare npedia pamdb phenolexplorer phytohub procardb respect sancdb streptomedb swmd symmap tmdb tmmc tppt triforc unpd #metabolights plantcyc 

DATABASES_RECONVERT = biofacquim-reconvert metabolights-reconvert streptomedb-reconvert

DATABASES_REINTEGRATE = plantcyc-reintegrate triforc-reintegrate unpd-reintegrate

DATABASES_RESCRAPE = alkamid-rescrape biophytmol-rescrape carotenoiddb-rescrape knapsack-rescrape metabolights-rescrape mitishamba-rescrape nanpdb-rescrape npedia-rescrape phytohub-rescrape procardb-rescrape sancdb-rescrape swmd-rescrape tmdb-rescrape

.PHONY: ${DATABASES} ${DATABATES_RECONVERT} ${DATABASES_REINTEGRATE} ${DATABASES_RESCRAPE}

all: ${DATABASES}
	@echo Cur dir: ${CURDIR}
	@echo PWD: ${PWD}
	@echo DB source path: ${DB_SOURCE_PATH}

databases: afrotryp alkamid biofacquim biophytmol carotenoiddb cmaup coconut cyanometdb dnp drduke etcm foodb inflamnat knapsack metabolights mibig mitishamba nanpdb npass npatlas npcare npedia pamdb phenolexplorer phytohub plantcyc procardb respect sancdb streptomedb swmd symmap tmdb tmmc tppt triforc unpd

databases-reconvert: ${DATABASES_RECONVERT}

databases-reintegrate: ${DATABASES_REINTEGRATE}

databases-rescrape: ${DATABASES_RESCRAPE}

${DATABASES}:%:${INTERIM_PATH}/%.tsv.zip

${INTERIM_PATH}/afrotryp.tsv.zip: ${DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip ${SRC_GATHERING_PATH}/afrotryp/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/afrotryp/standardizing.R

${INTERIM_PATH}/alkamid.tsv.zip: ${DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.zip ${DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.zip ${SRC_GATHERING_PATH}/alkamid/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/alkamid/standardizing.R

alkamid-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/alkamid/scraping.R

${INTERIM_PATH}/biofacquim.tsv.zip: ${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip ${SRC_GATHERING_PATH}/biofacquim/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/biofacquim/standardizing.R

biofacquim-reconvert: ${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip

${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip: ${DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf ${SRC_GATHERING_PATH}/biofacquim/converting.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/biofacquim/converting.R

${INTERIM_PATH}/biophytmol.tsv.zip: ${DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.zip ${SRC_GATHERING_PATH}/biophytmol/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/biophytmol/standardizing.R

biophytmol-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/biophytmol/scraping.R

${INTERIM_PATH}/carotenoiddb.tsv.zip: ${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip ${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R

carotenoiddb-rescrape: ${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip

${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip : ${DB_SOURCE_PATH}/carotenoiddb/Carotenoids_InChI_InChIKey.tsv ${SRC_GATHERING_PATH}/carotenoiddb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/carotenoiddb/scraping.R

${INTERIM_PATH}/cmaup.tsv.zip: ${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Ingredients_All.txt ${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plants.txt ${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt ${SRC_GATHERING_PATH}/cmaup/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/cmaup/standardizing.R

# maybe split both indeed else long
${INTERIM_PATH}/coconut.tsv.zip: ${DB_SOURCE_PATH}/coconut/COCONUT.sdf.zip ${DB_SOURCE_PATH}/coconut/coconutConverted.tsv.zip ${SRC_GATHERING_PATH}/coconut/converting.py ${SRC_GATHERING_PATH}/coconut/standardizing.R
	cd ${PWD}/src &&	python ${SRC_GATHERING_PATH}/coconut/converting.py  && Rscript ${SRC_GATHERING_PATH}/coconut/standardizing.R

${INTERIM_PATH}/cyanometdb.tsv.zip: ${DB_SOURCE_PATH}/cyanometdb/media-1.csv ${SRC_GATHERING_PATH}/cyanometdb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/cyanometdb/standardizing.R

${INTERIM_PATH}/dnp.tsv.zip: ${DB_SOURCE_PATH}/dnp/28_2/full_set.csv ${SRC_GATHERING_PATH}/dnp/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/dnp/standardizing.R

${INTERIM_PATH}/drduke.tsv.zip: ${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv ${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv ${SRC_GATHERING_PATH}/drduke/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/drduke/standardizing.R

${INTERIM_PATH}/etcm.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/etcm/data/*.csv) ${SRC_GATHERING_PATH}/etcm/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/etcm/standardizing.R

${INTERIM_PATH}/foodb.tsv.zip: ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Compound_copy.csv ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Content.csv ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Flavor.csv ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv ${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Reference.csv ${SRC_GATHERING_PATH}/foodb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/foodb/standardizing.R

${INTERIM_PATH}/inflamnat.tsv.zip: ${DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx ${SRC_GATHERING_PATH}/inflamnat/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/inflamnat/standardizing.R

${INTERIM_PATH}/knapsack.tsv.zip: ${DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.zip ${SRC_GATHERING_PATH}/knapsack/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/knapsack/standardizing.R

knapsack-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/knapsack/scraping.R

# PROCESS DONE FOR METABOLIGHTS IS HORROR... my apologies I can't write it another way
## to be changed

${INTERIM_PATH}/metabolights.tsv.zip: ${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip ${DB_SOURCE_PATH}/metabolights/metabolightsStudiesScraped.tsv.zip ${SRC_GATHERING_PATH}/metabolights/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/metabolights/standardizing.R

metabolights-reconvert: ${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip ${DB_SOURCE_PATH}/metabolights/metabolightsStudiesScraped.tsv.zip

${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/metabolights/studiesScraped/*.json) ${DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml ${SRC_GATHERING_PATH}/metabolights/prestandardizing.R ${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/metabolights/prestandardizing.R &&	Rscript ${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R

metabolights-rescrape: $(wildcard ${DB_SOURCE_PATH}/metabolights/studiesScraped/*.json)

$(wildcard ${DB_SOURCE_PATH}/metabolights/studiesScraped/*.json) : ${DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_studies.xml ${SRC_GATHERING_PATH}/metabolights/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/metabolights/scraping.R

${INTERIM_PATH}/mibig.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/mibig/mibig_json_2.0/*.json) ${SRC_GATHERING_PATH}/mibig/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/mibig/standardizing.R

${INTERIM_PATH}/mitishamba.tsv.zip: ${DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.zip ${SRC_GATHERING_PATH}/mitishamba/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/mitishamba/standardizing.R

mitishamba-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/mitishamba/scraping.R

${INTERIM_PATH}/nanpdb.tsv.zip: ${DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.zip ${SRC_GATHERING_PATH}/nanpdb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/nanpdb/standardizing.R

nanpdb-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/nanpdb/scraping.R

${INTERIM_PATH}/npass.tsv.zip: ${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt ${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_properties.txt ${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt ${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_species_pair.txt ${SRC_GATHERING_PATH}/npass/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/npass/standardizing.R

${INTERIM_PATH}/npatlas.tsv.zip: ${DB_SOURCE_PATH}/npatlas/np_atlas_2019_12.tsv ${SRC_GATHERING_PATH}/npatlas/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/npatlas/standardizing.R

${INTERIM_PATH}/npcare.tsv.zip: ${DB_SOURCE_PATH}/npcare/npcare.zip ${SRC_GATHERING_PATH}/npcare/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/npcare/standardizing.R

${INTERIM_PATH}/npedia.tsv.zip: ${DB_SOURCE_PATH}/npedia/npediaScraped.tsv.zip ${SRC_GATHERING_PATH}/npedia/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/npedia/standardizing.R

npedia-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/npedia/scraping.R

${INTERIM_PATH}/pamdb.tsv.zip: ${DB_SOURCE_PATH}/pamdb/PaMet.xlsx ${SRC_GATHERING_PATH}/pamdb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/pamdb/standardizing.R

${INTERIM_PATH}/phenolexplorer.tsv.zip: ${DB_SOURCE_PATH}/phenolexplorer/compounds-classification.csv ${DB_SOURCE_PATH}/phenolexplorer/compounds-structures.csv ${DB_SOURCE_PATH}/phenolexplorer/compounds.csv ${DB_SOURCE_PATH}/phenolexplorer/foods-classification.csv ${DB_SOURCE_PATH}/phenolexplorer/foods.csv ${DB_SOURCE_PATH}/phenolexplorer/metabolites-structures.csv ${DB_SOURCE_PATH}/phenolexplorer/metabolites.csv ${DB_SOURCE_PATH}/phenolexplorer/publications.csv ${DB_SOURCE_PATH}/phenolexplorer/composition-data.xlsx ${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R

${INTERIM_PATH}/phytohub.tsv.zip: ${DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.zip ${SRC_GATHERING_PATH}/phytohub/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/phytohub/standardizing.R

phytohub-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/phytohub/scraping.R

${INTERIM_PATH}/plantcyc.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/plantcyc/*.tsv.zip) ${SRC_GATHERING_PATH}/plantcyc/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/plantcyc/standardizing.R

# do not know how to proceed here... might be wrong
plantcyc-reintegrate: $(wildcard ${DB_SOURCE_PATH}/plantcyc/*.tsv.zip)

$(wildcard ${DB_SOURCE_PATH}/plantcyc/*.tsv.zip) : $(wildcard ${DB_SOURCE_PATH}/plantcyc/0_data/*/*/data/compounds.dat) ${SRC_GATHERING_PATH}/plantcyc/integrating.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/plantcyc/integrating.R

${INTERIM_PATH}/procardb.tsv.zip: ${DB_SOURCE_PATH}/procardb/procardbScraped.tsv.zip ${SRC_GATHERING_PATH}/procardb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/procardb/standardizing.R

prodardb-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/procardb/scraping.R

${INTERIM_PATH}/respect.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/respect/*.txt) ${SRC_GATHERING_PATH}/respect/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/respect/standardizing.R

${INTERIM_PATH}/sancdb.tsv.zip: ${DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.zip ${SRC_GATHERING_PATH}/sancdb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/sancdb/standardizing.R

sancdb-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/sancdb/scraping.R

${INTERIM_PATH}/streptomedb.tsv.zip: ${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip ${SRC_GATHERING_PATH}/streptomedb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/streptomedb/standardizing.R

streptomedb-reconvert: ${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip

${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip: ${DB_SOURCE_PATH}/streptomedb/streptomedb.sdf ${SRC_GATHERING_PATH}/streptomedb/converting.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/streptomedb/converting.R

${INTERIM_PATH}/swmd.tsv.zip: ${DB_SOURCE_PATH}/swmd/swmdScraped.tsv.zip ${SRC_GATHERING_PATH}/swmd/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/swmd/standardizing.R

swmd-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/swmd/scraping.R

${INTERIM_PATH}/symmap.tsv.zip: $(wildcard ${DB_SOURCE_PATH}/symmap/data/*.csv) ${SRC_GATHERING_PATH}/symmap/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/symmap/standardizing.R

${INTERIM_PATH}/tmdb.tsv.zip: ${DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.zip ${SRC_GATHERING_PATH}/tmdb/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/tmdb/standardizing.R

tmdb-rescrape:
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/tmdb/scraping.R

${INTERIM_PATH}/tmmc.tsv.zip: ${DB_SOURCE_PATH}/tmmc/compound.xlsx ${SRC_GATHERING_PATH}/tmmc/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/tmmc/standardizing.R

${INTERIM_PATH}/tppt.tsv.zip: ${DB_SOURCE_PATH}/tppt/TPPT_database.xlsx ${SRC_GATHERING_PATH}/tppt/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/tppt/standardizing.R

${INTERIM_PATH}/triforc.tsv.zip: ${DB_SOURCE_PATH}/triforc/triforcBis.tsv ${SRC_GATHERING_PATH}/triforc/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/triforc/standardizing.R

triforc-reintegrate: ${DB_SOURCE_PATH}/triforc/triforcBis.tsv

${DB_SOURCE_PATH}/triforc/triforcBis.tsv : ${DB_SOURCE_PATH}/triforc/triforcOriginal.tsv ${DB_SOURCE_PATH}/triforc/triforcToGet.tsv ${SRC_GATHERING_PATH}/triforc/prestandardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/triforc/prestandardizing.R

${INTERIM_PATH}/unpd.tsv.zip: ${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip ${SRC_GATHERING_PATH}/unpd/standardizing.R
	cd ${PWD}/src &&	Rscript ${SRC_GATHERING_PATH}/unpd/standardizing.R

unpd-reintegrate: ${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip

${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip: ${DB_SOURCE_PATH}/unpd/unpd_final.csv.zip ${DB_SOURCE_PATH}/unpd/UNPD_DB.csv.zip ${SRC_GATHERING_PATH}/unpd/integrating.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_PATH}/unpd/integrating.R
