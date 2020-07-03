include	../../../paths.mk

DATABASES = afrotryp	alkamid	biofacquim	biophytmol	carotenoiddb	cmaup	coconut	cyanometdb	datawarrior	dnp	drduke	etcm	foodb	inflamnat	knapsack	metabolights	mibig	mitishamba	nanpdb	npass	npatlas	npcare	npedia	nubbe	pamdb	pharmdbk	phenolexplorer	phytohub	procardb	plantcyc	respect	sancdb	streptomedb	swmd	symmap	tmdb	tmmc	tppt	triforc	unpd	wakankensaku

DATABASES_RECONVERT = biofacquim-reconvert	coconut-reconvert	metabolights-reconvert	streptomedb-reconvert

DATABASES_REINTEGRATE = plantcyc-reintegrate	triforc-reintegrate	unpd-reintegrate

DATABASES_RESCRAPE = alkamid-rescrape	biophytmol-rescrape	carotenoiddb-rescrape	knapsack-rescrape	metabolights-rescrape	mitishamba-rescrape	nanpdb-rescrape	npedia-rescrape	phytohub-rescrape	procardb-rescrape	sancdb-rescrape	swmd-rescrape	tmdb-rescrape	wakankensaku-rescrape	# nubbe-rescrape pharmdbk-rescrape not done

.PHONY:	${DATABASES}	${DATABATES_RECONVERT}	${DATABASES_REINTEGRATE}	${DATABASES_RESCRAPE}

all: ${DATABASES}
	@echo Cur dir: ${CURDIR}
	@echo PWD: ${PWD}
	@echo DB source path: ${DB_SOURCE_PATH}

databases:	${DATABASES}

databases-reconvert:	${DATABASES_RECONVERT}

databases-reintegrate:	${DATABASES_REINTEGRATE}

databases-rescrape:	${DATABASES_RESCRAPE}

${DATABASES}:%:${INTERIM_DB_PATH}/%.tsv.zip

afrotryp:	${INTERIM_DB_PATH}/afrotryp.tsv.zip
${INTERIM_DB_PATH}/afrotryp.tsv.zip: ${DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip	${SRC_GATHERING_PATH}/afrotryp/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/afrotryp/standardizing.R

alkamid:	${INTERIM_DB_PATH}/alkamid.tsv.zip
${INTERIM_DB_PATH}/alkamid.tsv.zip:	${DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.zip	${DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.zip	${SRC_GATHERING_PATH}/alkamid/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/alkamid/standardizing.R

alkamid-rescrape:	${DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.zip	${DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.zip
${DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.zip	${DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/alkamid/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/alkamid/scraping.R

biofacquim:	${INTERIM_DB_PATH}/biofacquim.tsv.zip
${INTERIM_DB_PATH}/biofacquim.tsv.zip:	${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip	${SRC_GATHERING_PATH}/biofacquim/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/biofacquim/standardizing.R

biofacquim-reconvert:	${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip
${DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.zip:	${DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf	${SRC_GATHERING_PATH}/biofacquim/converting.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/biofacquim/converting.R

biophytmol:	${INTERIM_DB_PATH}/biophytmol.tsv.zip
${INTERIM_DB_PATH}/biophytmol.tsv.zip:	${DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.zip	${SRC_GATHERING_PATH}/biophytmol/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/biophytmol/standardizing.R

biophytmol-rescrape:	${DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.zip
${DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/biophytmol/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/biophytmol/scraping.R

carotenoiddb:	${INTERIM_DB_PATH}/carotenoiddb.tsv.zip
${INTERIM_DB_PATH}/carotenoiddb.tsv.zip:	${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip	${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/carotenoiddb/standardizing.R

carotenoiddb-rescrape:	${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip
${DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.zip:	${DB_SOURCE_PATH}/carotenoiddb/Carotenoids_InChI_InChIKey.tsv	# ${SRC_GATHERING_PATH}/carotenoiddb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/carotenoiddb/scraping.R

cmaup:	${INTERIM_DB_PATH}/cmaup.tsv.zip
${INTERIM_DB_PATH}/cmaup.tsv.zip:	${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Ingredients_All.txt	${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plants.txt	${DB_SOURCE_PATH}/cmaup/CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt	${SRC_GATHERING_PATH}/cmaup/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/cmaup/standardizing.R

coconut:	${INTERIM_DB_PATH}/coconut.tsv.zip
${INTERIM_DB_PATH}/coconut.tsv.zip:	${DB_SOURCE_PATH}/coconut/coconutConverted.tsv.zip	${SRC_GATHERING_PATH}/coconut/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/coconut/standardizing.R

coconut-reconvert:	${DB_SOURCE_PATH}/coconut/coconutConverted.tsv.zip
${DB_SOURCE_PATH}/coconut/coconutConverted.tsv.zip:	${DB_SOURCE_PATH}/coconut/COCONUT_DB.sdf	${SRC_GATHERING_PATH}/coconut/converting.py
	cd	${PWD}/src	&&	python	${SRC_GATHERING_PATH}/coconut/converting.py

cyanometdb:	${INTERIM_DB_PATH}/cyanometdb.tsv.zip
${INTERIM_DB_PATH}/cyanometdb.tsv.zip:	${DB_SOURCE_PATH}/cyanometdb/media-1.csv	${SRC_GATHERING_PATH}/cyanometdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/cyanometdb/standardizing.R

datawarrior:	${INTERIM_DB_PATH}/datawarrior.tsv.zip
${INTERIM_DB_PATH}/datawarrior.tsv.zip:	${DB_SOURCE_PATH}/datawarrior/NaturalProducts.txt	${SRC_GATHERING_PATH}/datawarrior/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/datawarrior/standardizing.R

dnp:	${INTERIM_DB_PATH}/dnp.tsv.zip
${INTERIM_DB_PATH}/dnp.tsv.zip:	${DB_SOURCE_PATH}/dnp/28_2/full_set.csv	${SRC_GATHERING_PATH}/dnp/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/dnp/standardizing.R

drduke:	${INTERIM_DB_PATH}/drduke.tsv.zip
${INTERIM_DB_PATH}/drduke.tsv.zip:	${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv	${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv	${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv	${DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv	${SRC_GATHERING_PATH}/drduke/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/drduke/standardizing.R

etcm:	${INTERIM_DB_PATH}/etcm.tsv.zip
${INTERIM_DB_PATH}/etcm.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/etcm/data/*.csv)	${SRC_GATHERING_PATH}/etcm/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/etcm/standardizing.R

foodb:	${INTERIM_DB_PATH}/foodb.tsv.zip
${INTERIM_DB_PATH}/foodb.tsv.zip:	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Compound_copy.csv	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Content.csv	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Flavor.csv	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv	${DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Reference.csv	${SRC_GATHERING_PATH}/foodb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/foodb/standardizing.R

inflamnat:	${INTERIM_DB_PATH}/inflamnat.tsv.zip
${INTERIM_DB_PATH}/inflamnat.tsv.zip:	${DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx	${SRC_GATHERING_PATH}/inflamnat/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/inflamnat/standardizing.R

knapsack:	${INTERIM_DB_PATH}/knapsack.tsv.zip
${INTERIM_DB_PATH}/knapsack.tsv.zip:	${DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.zip	${SRC_GATHERING_PATH}/knapsack/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/knapsack/standardizing.R

knapsack-rescrape:	${DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.zip
${DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/knapsack/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/knapsack/scraping.R

# PROCESS DONE FOR METABOLIGHTS IS HORROR... my apologies I can't write it another way
## to be changed

metabolights:	${INTERIM_DB_PATH}/metabolights.tsv.zip
${INTERIM_DB_PATH}/metabolights.tsv.zip:	${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip	${DB_SOURCE_PATH}/metabolights/metabolightsStudiesScraped.tsv.zip	${SRC_GATHERING_PATH}/metabolights/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/metabolights/standardizing.R

metabolights-reconvert:	${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip
${DB_SOURCE_PATH}/metabolights/metabolightsPrecleaned.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/metabolights/studiesScraped/*.json)	${DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml	# ${SRC_GATHERING_PATH}/metabolights/prestandardizing.R	${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/metabolights/prestandardizing.R	&&	Rscript	${SRC_GATHERING_PATH}/metabolights/standardizingStudies.R

metabolights-rescrape:	${DB_SOURCE_PATH}/metabolights/metabolightsStudiesScraped.tsv.zip
${DB_SOURCE_PATH}/metabolights/metabolightsStudiesScraped.tsv.zip:	${DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_studies.xml	# ${SRC_GATHERING_PATH}/metabolights/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/metabolights/scraping.R

${INTERIM_DB_PATH}/mibig.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/mibig/mibig_json_2.0/*.json)	${SRC_GATHERING_PATH}/mibig/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/mibig/standardizing.R

mitishamba:	${INTERIM_DB_PATH}/mitishamba.tsv.zip
${INTERIM_DB_PATH}/mitishamba.tsv.zip:	${DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.zip	${SRC_GATHERING_PATH}/mitishamba/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/mitishamba/standardizing.R

mitishamba-rescrape:	${DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.zip
${DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/mitishamba/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/mitishamba/scraping.R

nanpdb:	${INTERIM_DB_PATH}/nanpdb.tsv.zip
${INTERIM_DB_PATH}/nanpdb.tsv.zip:	${DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.zip	${SRC_GATHERING_PATH}/nanpdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/nanpdb/standardizing.R

nanpdb-rescrape:	${DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.zip
${DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/nanpdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/nanpdb/scraping.R

npass:	${INTERIM_DB_PATH}/npass.tsv.zip
${INTERIM_DB_PATH}/npass.tsv.zip:	${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt	${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_properties.txt	${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt	${DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_species_pair.txt	${SRC_GATHERING_PATH}/npass/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/npass/standardizing.R

npatlas:	${INTERIM_DB_PATH}/npatlas.tsv.zip
${INTERIM_DB_PATH}/npatlas.tsv.zip:	${DB_SOURCE_PATH}/npatlas/np_atlas_2019_12.tsv	${SRC_GATHERING_PATH}/npatlas/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/npatlas/standardizing.R

npcare:	${INTERIM_DB_PATH}/npcare.tsv.zip
${INTERIM_DB_PATH}/npcare.tsv.zip:	${DB_SOURCE_PATH}/npcare/npcare.zip	${SRC_GATHERING_PATH}/npcare/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/npcare/standardizing.R

npedia:	${INTERIM_DB_PATH}/npedia.tsv.zip
${INTERIM_DB_PATH}/npedia.tsv.zip:	${DB_SOURCE_PATH}/npedia/npediaScraped.tsv.zip	${SRC_GATHERING_PATH}/npedia/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/npedia/standardizing.R

npedia-rescrape:	${DB_SOURCE_PATH}/npedia/npediaScraped.tsv.zip
${DB_SOURCE_PATH}/npedia/npediaScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/npedia/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/npedia/scraping.R

nubbe:	${INTERIM_DB_PATH}/nubbe.tsv.zip
${INTERIM_DB_PATH}/nubbe.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/nubbe/nubbe_raw/*.xml)	${SRC_GATHERING_PATH}/nubbe/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/nubbe/standardizing.R

# nubbe-rescrape not done

pamdb:	${INTERIM_DB_PATH}/pamdb.tsv.zip
${INTERIM_DB_PATH}/pamdb.tsv.zip:	${DB_SOURCE_PATH}/pamdb/PaMet.xlsx	${SRC_GATHERING_PATH}/pamdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/pamdb/standardizing.R

pharmdbk:	${INTERIM_DB_PATH}/pharmdbk.tsv.zip
${INTERIM_DB_PATH}/pharmdbk.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/pharmdbk/pharmdbk_raw/*.json.gz)	${SRC_GATHERING_PATH}/pharmdbk/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/pharmdbk/standardizing.R

# pharmdbk-rescrape not done

phenolexplorer:	${INTERIM_DB_PATH}/phenolexplorer.tsv.zip
${INTERIM_DB_PATH}/phenolexplorer.tsv.zip:	${DB_SOURCE_PATH}/phenolexplorer/compounds-classification.csv	${DB_SOURCE_PATH}/phenolexplorer/compounds-structures.csv	${DB_SOURCE_PATH}/phenolexplorer/compounds.csv	${DB_SOURCE_PATH}/phenolexplorer/foods-classification.csv	${DB_SOURCE_PATH}/phenolexplorer/foods.csv	${DB_SOURCE_PATH}/phenolexplorer/metabolites-structures.csv	${DB_SOURCE_PATH}/phenolexplorer/metabolites.csv	${DB_SOURCE_PATH}/phenolexplorer/publications.csv	${DB_SOURCE_PATH}/phenolexplorer/composition-data.xlsx	${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/phenolexplorer/standardizing.R

phytohub:	${INTERIM_DB_PATH}/phytohub.tsv.zip
${INTERIM_DB_PATH}/phytohub.tsv.zip:	${DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.zip	${SRC_GATHERING_PATH}/phytohub/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/phytohub/standardizing.R

phytohub-rescrape:	${DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.zip
${DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/phytohub/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/phytohub/scraping.R

plantcyc:	${INTERIM_DB_PATH}/plantcyc.tsv.zip
${INTERIM_DB_PATH}/plantcyc.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/plantcyc/*.tsv.zip)	${SRC_GATHERING_PATH}/plantcyc/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/plantcyc/standardizing.R

plantcyc-reintegrate:	$(wildcard	${DB_SOURCE_PATH}/plantcyc/*.tsv.zip)	
$(wildcard	${DB_SOURCE_PATH}/plantcyc/*.tsv.zip):	$(wildcard	${DB_SOURCE_PATH}/plantcyc/0_data/*/*/data/compounds.dat)
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/plantcyc/integrating.R

procardb:	${INTERIM_DB_PATH}/procardb.tsv.zip
${INTERIM_DB_PATH}/procardb.tsv.zip:	${DB_SOURCE_PATH}/procardb/procardbScraped.tsv.zip	${SRC_GATHERING_PATH}/procardb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/procardb/standardizing.R

prodardb-rescrape:	${DB_SOURCE_PATH}/procardb/procardbScraped.tsv.zip
${DB_SOURCE_PATH}/procardb/procardbScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/procardb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/procardb/scraping.R

respect:	${INTERIM_DB_PATH}/respect.tsv.zip
${INTERIM_DB_PATH}/respect.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/respect/*.txt)	${SRC_GATHERING_PATH}/respect/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/respect/standardizing.R

sancdb:	${INTERIM_DB_PATH}/sancdb.tsv.zip
${INTERIM_DB_PATH}/sancdb.tsv.zip:	${DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.zip	${SRC_GATHERING_PATH}/sancdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/sancdb/standardizing.R

sancdb-rescrape:${DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.zip
${DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/sancdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/sancdb/scraping.R

streptomedb:	${INTERIM_DB_PATH}/streptomedb.tsv.zip
${INTERIM_DB_PATH}/streptomedb.tsv.zip:	${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip	${SRC_GATHERING_PATH}/streptomedb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/streptomedb/standardizing.R

streptomedb-reconvert:	${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip
${DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.zip:	${DB_SOURCE_PATH}/streptomedb/streptomedb.sdf	${SRC_GATHERING_PATH}/streptomedb/converting.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/streptomedb/converting.R

swmd:	${INTERIM_DB_PATH}/swmd.tsv.zip
${INTERIM_DB_PATH}/swmd.tsv.zip:	${DB_SOURCE_PATH}/swmd/swmdScraped.tsv.zip	${SRC_GATHERING_PATH}/swmd/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/swmd/standardizing.R

swmd-rescrape:	${DB_SOURCE_PATH}/swmd/swmdScraped.tsv.zip
${DB_SOURCE_PATH}/swmd/swmdScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/swmd/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/swmd/scraping.R

symmap:	${INTERIM_DB_PATH}/symmap.tsv.zip
${INTERIM_DB_PATH}/symmap.tsv.zip:	$(wildcard	${DB_SOURCE_PATH}/symmap/data/*.csv)	${SRC_GATHERING_PATH}/symmap/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/symmap/standardizing.R

tmdb:	${INTERIM_DB_PATH}/tmdb.tsv.zip
${INTERIM_DB_PATH}/tmdb.tsv.zip:	${DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.zip	${SRC_GATHERING_PATH}/tmdb/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/tmdb/standardizing.R

tmdb-rescrape:	${DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.zip
${DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/tmdb/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/tmdb/scraping.R

tmmc:	${INTERIM_DB_PATH}/tmmc.tsv.zip
${INTERIM_DB_PATH}/tmmc.tsv.zip:	${DB_SOURCE_PATH}/tmmc/compound.xlsx	${SRC_GATHERING_PATH}/tmmc/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/tmmc/standardizing.R

tppt:	${INTERIM_DB_PATH}/tppt.tsv.zip
${INTERIM_DB_PATH}/tppt.tsv.zip:	${DB_SOURCE_PATH}/tppt/TPPT_database.xlsx	${SRC_GATHERING_PATH}/tppt/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/tppt/standardizing.R

triforc:	${INTERIM_DB_PATH}/triforc.tsv.zip
${INTERIM_DB_PATH}/triforc.tsv.zip:	${DB_SOURCE_PATH}/triforc/triforcBis.tsv	${SRC_GATHERING_PATH}/triforc/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/triforc/standardizing.R

triforc-reintegrate:	${DB_SOURCE_PATH}/triforc/triforcBis.tsv
${DB_SOURCE_PATH}/triforc/triforcBis.tsv:
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/triforc/prestandardizing.R

unpd:	${INTERIM_DB_PATH}/unpd.tsv.zip
${INTERIM_DB_PATH}/unpd.tsv.zip:	${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip	${SRC_GATHERING_PATH}/unpd/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/unpd/standardizing.R

unpd-reintegrate:	${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip
${DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.zip:
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/unpd/integrating.R

wakankensaku:	${INTERIM_DB_PATH}/wakankensaku.tsv.zip
${INTERIM_DB_PATH}/wakankensaku.tsv.zip:	${DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.zip	${SRC_GATHERING_PATH}/wakankensaku/standardizing.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/wakankensaku/standardizing.R

wakankensaku-rescrape:	${DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.zip
${DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.zip:	# ${SRC_GATHERING_PATH}/wakankensaku/scraping.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_PATH}/wakankensaku/scraping.R
