include	../../../paths.mk

DATABASES = afrotryp alkamid biofacquim biophytmol carotenoiddb coconut cyanometdb datawarrior dnp drduke foodb inflamnat knapsack metabolights mibig mitishamba nanpdb napralert npass npatlas npcare npedia nubbe pamdb phenolexplorer phytohub procardb respect sancdb streptomedb swmd tmdb tmmc tppt unpd wakankensaku # cmaup etcm pharmdbk plantcyc symmap triforc

DATABASES_RECONVERT = biofacquim-reconvert coconut-reconvert streptomedb-reconvert

DATABASES_REINTEGRATE = unpd-reintegrate # plantcyc-reintegrate 

DATABASES_RESCRAPE = alkamid-rescrape biophytmol-rescrape carotenoiddb-rescrape knapsack-rescrape mitishamba-rescrape nanpdb-rescrape npedia-rescrape phytohub-rescrape procardb-rescrape sancdb-rescrape swmd-rescrape tmdb-rescrape wakankensaku-rescrape # nubbe-rescrape pharmdbk-rescrape not done

.PHONY: ${DATABASES} ${DATABATES_RECONVERT} ${DATABASES_REINTEGRATE} ${DATABASES_RESCRAPE}
.PHONY: gathering-databases gathering-databases-rescrape gathering-databases-reintegrate gathering-databases-reconvert
.PRECIOUS: %.tsv %.zip %.json %.gz

all: ${DATABASES}
	@echo Cur dir: ${CURDIR}
	@echo PWD: ${PWD}
	@echo DB source path: ${EXTERNAL_DB_SOURCE_PATH}

gathering-databases: ${DATABASES}

gathering-databases-reconvert: ${DATABASES_RECONVERT}

gathering-databases-reintegrate: ${DATABASES_REINTEGRATE}

gathering-databases-rescrape: ${DATABASES_RESCRAPE}

${DATABASES}:%:${INTERIM_DB_PATH}/%.tsv.gz

afrotryp: ${INTERIM_DB_PATH}/afrotryp.tsv.gz
${INTERIM_DB_PATH}/afrotryp.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip ${SRC_GATHERING_DB_PATH}/afrotryp/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/afrotryp/standardizing.R

alkamid: ${INTERIM_DB_PATH}/alkamid.tsv.gz
${INTERIM_DB_PATH}/alkamid.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/alkamid/*.tsv.gz) ${SRC_GATHERING_DB_PATH}/alkamid/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/alkamid/standardizing.R

alkamid-rescrape: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/alkamid/*.tsv.gz)
$(wildcard ${EXTERNAL_DB_SOURCE_PATH}/alkamid/*.tsv.gz): # ${SRC_GATHERING_DB_PATH}/alkamid/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/alkamid/scraping.R

biofacquim: ${INTERIM_DB_PATH}/biofacquim.tsv.gz
${INTERIM_DB_PATH}/biofacquim.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz ${SRC_GATHERING_DB_PATH}/biofacquim/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/biofacquim/standardizing.R

biofacquim-reconvert: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf ${SRC_GATHERING_DB_PATH}/biofacquim/converting.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/biofacquim/converting.R

biophytmol: ${INTERIM_DB_PATH}/biophytmol.tsv.gz
${INTERIM_DB_PATH}/biophytmol.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/biophytmol/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/biophytmol/standardizing.R

biophytmol-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/biophytmol/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/biophytmol/scraping.R

carotenoiddb: ${INTERIM_DB_PATH}/carotenoiddb.tsv.gz
${INTERIM_DB_PATH}/carotenoiddb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/carotenoiddb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/carotenoiddb/standardizing.R

carotenoiddb-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/Carotenoids_InChI_InChIKey.tsv ${SRC_GATHERING_DB_PATH}/carotenoiddb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/carotenoiddb/scraping.R

# cmaup: ${INTERIM_DB_PATH}/cmaup.tsv.gz
# ${INTERIM_DB_PATH}/cmaup.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Ingredients_All.txt ${EXTERNAL_DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plants.txt ${EXTERNAL_DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt ${SRC_GATHERING_DB_PATH}/cmaup/standardizing.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/cmaup/standardizing.R

coconut: ${INTERIM_DB_PATH}/coconut.tsv.gz
${INTERIM_DB_PATH}/coconut.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz ${SRC_GATHERING_DB_PATH}/coconut/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/coconut/standardizing.R

coconut-reconvert: ${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/coconut/COCONUT_DB.sdf ${SRC_GATHERING_DB_PATH}/coconut/converting.py
	cd	${PWD}/src	&&	python	${SRC_GATHERING_DB_PATH}/coconut/converting.py

cyanometdb: ${INTERIM_DB_PATH}/cyanometdb.tsv.gz
${INTERIM_DB_PATH}/cyanometdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/cyanometdb/media-1.csv ${SRC_GATHERING_DB_PATH}/cyanometdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/cyanometdb/standardizing.R

datawarrior: ${INTERIM_DB_PATH}/datawarrior.tsv.gz
${INTERIM_DB_PATH}/datawarrior.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/datawarrior/NaturalProducts.txt ${SRC_GATHERING_DB_PATH}/datawarrior/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/datawarrior/standardizing.R

dnp: ${INTERIM_DB_PATH}/dnp.tsv.gz
${INTERIM_DB_PATH}/dnp.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/dnp/29_2/full_set.csv ${SRC_GATHERING_DB_PATH}/dnp/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/dnp/standardizing.R

drduke: ${INTERIM_DB_PATH}/drduke.tsv.gz
${INTERIM_DB_PATH}/drduke.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv ${SRC_GATHERING_DB_PATH}/drduke/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/drduke/standardizing.R

# etcm: ${INTERIM_DB_PATH}/etcm.tsv.gz
# ${INTERIM_DB_PATH}/etcm.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/etcm/data.zip ${SRC_GATHERING_DB_PATH}/etcm/standardizing.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/etcm/standardizing.R

foodb: ${INTERIM_DB_PATH}/foodb.tsv.gz
${INTERIM_DB_PATH}/foodb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Compound_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Content.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Flavor.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Reference.csv ${SRC_GATHERING_DB_PATH}/foodb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/foodb/standardizing.R

inflamnat: ${INTERIM_DB_PATH}/inflamnat.tsv.gz
${INTERIM_DB_PATH}/inflamnat.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx ${SRC_GATHERING_DB_PATH}/inflamnat/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/inflamnat/standardizing.R

knapsack: ${INTERIM_DB_PATH}/knapsack.tsv.gz
${INTERIM_DB_PATH}/knapsack.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/knapsack/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/knapsack/standardizing.R

knapsack-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/knapsack/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/knapsack/scraping.R

metabolights: ${INTERIM_DB_PATH}/metabolights.tsv.gz
${INTERIM_DB_PATH}/metabolights.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml ${SRC_GATHERING_DB_PATH}/metabolights/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/metabolights/standardizing.R

mibig: ${INTERIM_DB_PATH}/mibig.tsv.gz
${INTERIM_DB_PATH}/mibig.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/mibig/data.zip ${SRC_GATHERING_DB_PATH}/mibig/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/mibig/standardizing.R

mitishamba: ${INTERIM_DB_PATH}/mitishamba.tsv.gz
${INTERIM_DB_PATH}/mitishamba.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/mitishamba/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/mitishamba/standardizing.R

mitishamba-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/mitishamba/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/mitishamba/scraping.R

nanpdb: ${INTERIM_DB_PATH}/nanpdb.tsv.gz
${INTERIM_DB_PATH}/nanpdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/nanpdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/nanpdb/standardizing.R

nanpdb-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz:	# ${SRC_GATHERING_DB_PATH}/nanpdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/nanpdb/scraping.R

napralert: ${INTERIM_DB_PATH}/napralert.tsv.gz
${INTERIM_DB_PATH}/napralert.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert_matched_final_unified.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert_matched_final_unified.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert.tsv.gz ${SRC_GATHERING_DB_PATH}/napralert/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/napralert/standardizing.R

npass: ${INTERIM_DB_PATH}/npass.tsv.gz
${INTERIM_DB_PATH}/npass.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_properties.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_species_pair.txt ${SRC_GATHERING_DB_PATH}/npass/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/npass/standardizing.R

npatlas: ${INTERIM_DB_PATH}/npatlas.tsv.gz
${INTERIM_DB_PATH}/npatlas.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npatlas/NPAtlas_download.tsv ${SRC_GATHERING_DB_PATH}/npatlas/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/npatlas/standardizing.R

npcare: ${INTERIM_DB_PATH}/npcare.tsv.gz
${INTERIM_DB_PATH}/npcare.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npcare/npcare.zip ${SRC_GATHERING_DB_PATH}/npcare/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/npcare/standardizing.R

npedia: ${INTERIM_DB_PATH}/npedia.tsv.gz
${INTERIM_DB_PATH}/npedia.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/npedia/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/npedia/standardizing.R

npedia-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/npedia/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/npedia/scraping.R

nubbe: ${INTERIM_DB_PATH}/nubbe.tsv.gz
${INTERIM_DB_PATH}/nubbe.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/nubbe/nubbe_raw/*.xml) ${SRC_GATHERING_DB_PATH}/nubbe/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/nubbe/standardizing.R

# nubbe-rescrape not done

pamdb: ${INTERIM_DB_PATH}/pamdb.tsv.gz
${INTERIM_DB_PATH}/pamdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/pamdb/PaMet.xlsx ${SRC_GATHERING_DB_PATH}/pamdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/pamdb/standardizing.R

# pharmdbk: ${INTERIM_DB_PATH}/pharmdbk.tsv.gz
# ${INTERIM_DB_PATH}/pharmdbk.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/pharmdbk/pharmdbk_raw/*.json.gz) ${SRC_GATHERING_DB_PATH}/pharmdbk/standardizing.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/pharmdbk/standardizing.R

# pharmdbk-rescrape not done

phenolexplorer: ${INTERIM_DB_PATH}/phenolexplorer.tsv.gz
${INTERIM_DB_PATH}/phenolexplorer.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/compounds-classification.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/compounds-structures.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/compounds.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods-classification.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/metabolites-structures.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/metabolites.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/publications.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/composition-data.xlsx ${SRC_GATHERING_DB_PATH}/phenolexplorer/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/phenolexplorer/standardizing.R

phytohub: ${INTERIM_DB_PATH}/phytohub.tsv.gz
${INTERIM_DB_PATH}/phytohub.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/phytohub/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/phytohub/standardizing.R

phytohub-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/phytohub/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/phytohub/scraping.R

# plantcyc: ${INTERIM_DB_PATH}/plantcyc.tsv.gz
# ${INTERIM_DB_PATH}/plantcyc.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/plantcyc/*.tsv.gz) ${SRC_GATHERING_DB_PATH}/plantcyc/standardizing.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/plantcyc/standardizing.R

# plantcyc-reintegrate: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/plantcyc/*.tsv.gz)	
# $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/plantcyc/*.tsv.gz): # $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/plantcyc/0_data/*/*/data/compounds.dat) ${SRC_GATHERING_DB_PATH}/plantcyc/integrating.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/plantcyc/integrating.R

procardb: ${INTERIM_DB_PATH}/procardb.tsv.gz
${INTERIM_DB_PATH}/procardb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/procardb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/procardb/standardizing.R

procardb-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz:	# ${SRC_GATHERING_DB_PATH}/procardb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/procardb/scraping.R

respect: ${INTERIM_DB_PATH}/respect.tsv.gz
${INTERIM_DB_PATH}/respect.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/respect/*.txt) ${SRC_GATHERING_DB_PATH}/respect/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/respect/standardizing.R

sancdb: ${INTERIM_DB_PATH}/sancdb.tsv.gz
${INTERIM_DB_PATH}/sancdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/sancdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/sancdb/standardizing.R

sancdb-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/sancdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/sancdb/scraping.R

streptomedb: ${INTERIM_DB_PATH}/streptomedb.tsv.gz
${INTERIM_DB_PATH}/streptomedb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz ${SRC_GATHERING_DB_PATH}/streptomedb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/streptomedb/standardizing.R

streptomedb-reconvert: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.sdf ${SRC_GATHERING_DB_PATH}/streptomedb/converting.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/streptomedb/converting.R

swmd: ${INTERIM_DB_PATH}/swmd.tsv.gz
${INTERIM_DB_PATH}/swmd.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/swmd/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/swmd/standardizing.R

swmd-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/swmd/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/swmd/scraping.R

# symmap: ${INTERIM_DB_PATH}/symmap.tsv.gz
# ${INTERIM_DB_PATH}/symmap.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/symmap/data.zip ${SRC_GATHERING_DB_PATH}/symmap/standardizing.R
# 	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/symmap/standardizing.R

tmdb: ${INTERIM_DB_PATH}/tmdb.tsv.gz
${INTERIM_DB_PATH}/tmdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/tmdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/tmdb/standardizing.R

tmdb-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz:	# ${SRC_GATHERING_DB_PATH}/tmdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/tmdb/scraping.R

tmmc: ${INTERIM_DB_PATH}/tmmc.tsv.gz
${INTERIM_DB_PATH}/tmmc.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx ${SRC_GATHERING_DB_PATH}/tmmc/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/tmmc/standardizing.R

tppt: ${INTERIM_DB_PATH}/tppt.tsv.gz
${INTERIM_DB_PATH}/tppt.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tppt/TPPT_database.xlsx ${SRC_GATHERING_DB_PATH}/tppt/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/tppt/standardizing.R

# triforc: ${INTERIM_DB_PATH}/triforc.tsv.gz
# ${INTERIM_DB_PATH}/triforc.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/triforc/triforcBis.tsv ${SRC_GATHERING_DB_PATH}/triforc/standardizing.R
#	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/triforc/standardizing.R

unpd: ${INTERIM_DB_PATH}/unpd.tsv.gz
${INTERIM_DB_PATH}/unpd.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz ${SRC_GATHERING_DB_PATH}/unpd/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/unpd/standardizing.R

unpd-reintegrate: ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/unpd/UNPD_DB.csv.zip ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpd_final.csv.zip ${SRC_GATHERING_DB_PATH}/unpd/integrating.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/unpd/integrating.R

wakankensaku: ${INTERIM_DB_PATH}/wakankensaku.tsv.gz
${INTERIM_DB_PATH}/wakankensaku.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/wakankensaku/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/wakankensaku/standardizing.R

wakankensaku-rescrape: ${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/wakankensaku/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_DB_PATH}/wakankensaku/scraping.R
