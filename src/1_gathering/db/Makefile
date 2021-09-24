DATABASES_ACCESSIBLE = creating-directory-db afrotryp biofacquim coconut cyanometdb datawarrior dianatdb drduke metabolights mibig npass npatlas npcare pamdb phenolexplorer respect streptomedb tmmc tppt# cmaup etcm pharmdbk plantcyc symmap triforc

DATABASES_SEMI = foodb-download foodb-download-modified foodb inflamnat-download inflamnat napralert nubbe unpd-download unpd-integrate unpd

DATABASES_CLOSED = creating-directory-db antibase-convert antibase antimarin-convert antimarin dnp

DATABASES_DOWNLOAD = creating-directory-db afrotryp-download biofacquim-download coconut-download coconut-convert cyanometdb-download dianatdb-download drduke-download metabolights-download mibig-download napralert-download npass-download npatlas-download npcare-download pamdb-download phenolexplorer-download respect-download streptomedb-download streptomedb-convert tmmc-download tppt-download

DATABASES_DOWNLOAD_MODIFIED = afrotryp-download-modified datawarrior-download-modified

DATABASES_SCRAPE = creating-directory-db alkamid-scrape alkamid biophytmol-scrape biofacquim-convert biophytmol carotenoiddb-scrape carotenoiddb knapsack-scrape knapsack mitishamba-scrape mitishamba nanpdb-scrape nanpdb npedia-scrape npedia phytohub-scrape phytohub procardb-scrape procardb sancdb-scrape sancdb swmd-scrape swmd tmdb-scrape tmdb wakankensaku-scrape wakankensaku# nubbe-scrape not done

.PHONY: ${DATABASES} ${DATABASES_DOWNLOAD} ${DATABASES_DOWNLOAD_MODIFIED} ${DATABASES_SCRAPE}
.PHONY: gathering-databases-accessible gathering-databases-scrape gathering-databases-convert creating-directory-db
.PRECIOUS: %.tsv %.zip %.json %.gz %.doc %.sdf

all: ${DATABASES}
	@echo Cur dir: ${CURDIR}
	@echo PWD: ${PWD}
	@echo DB source path: ${EXTERNAL_DB_SOURCE_PATH}

gathering-databases-accessible: ${DATABASES_ACCESSIBLE}

gathering-databases-semi: ${DATABASES_SEMI}

gathering-databases-closed: ${DATABASES_CLOSED}

gathering-databases-download: ${DATABASES_DOWNLOAD}

gathering-databases-download-modified: ${DATABASES_DOWNLOAD_MODIFIED}

gathering-databases-scrape: ${DATABASES_SCRAPE}

${DATABASES}:%:${INTERIM_DB_PATH}/%.tsv.gz

creating-directory-db:
	mkdir -p ${INTERIM_DB_PATH}

afrotryp: ${INTERIM_DB_PATH}/afrotryp.tsv.gz
${INTERIM_DB_PATH}/afrotryp.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip ${SRC_GATHERING_DB_PATH}/afrotryp/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/afrotryp/standardizing.R

afrotryp-download: ${EXTERNAL_DB_SOURCE_PATH}/afrotryp/44_2016_1764_MOESM1_ESM.doc
${EXTERNAL_DB_SOURCE_PATH}/afrotryp/44_2016_1764_MOESM1_ESM.doc: ${SRC_GATHERING_DB_PATH}/afrotryp/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/afrotryp/downloading.sh

afrotryp-download-modified: ${EXTERNAL_DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip
${EXTERNAL_DB_SOURCE_PATH}/afrotryp/afrotryp.tsv.zip: ${SRC_GATHERING_DB_PATH}/afrotryp/downloading_modified.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/afrotryp/downloading_modified.sh

alkamid: ${INTERIM_DB_PATH}/alkamid.tsv.gz
${INTERIM_DB_PATH}/alkamid.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/alkamid/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/alkamid/standardizing.R

alkamid-scrape: # ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.gz
# ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidScraped.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/alkamid/alkamidRefScraped.tsv.gz: ${SRC_GATHERING_DB_PATH}/alkamid/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/alkamid/scraping.R

antibase: ${INTERIM_DB_PATH}/antibase.tsv.gz
${INTERIM_DB_PATH}/antibase.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.smi ${SRC_GATHERING_DB_PATH}/antibase/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/antibase/standardizing.R

antibase-convert: ${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.smi ${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.smi ${EXTERNAL_DB_SOURCE_PATH}/antibase/antibaseConverted.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/antibase/ANTIBASE_2012_FORM2.sdf ${SRC_GATHERING_DB_PATH}/antibase/converting.py
	cd ${PWD}/src && python ${SRC_GATHERING_DB_PATH}/antibase/converting.py

antimarin: ${INTERIM_DB_PATH}/antimarin.tsv.gz
${INTERIM_DB_PATH}/antimarin.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.smi ${SRC_GATHERING_DB_PATH}/antimarin/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/antimarin/standardizing.R

antimarin-convert: ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.smi ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.smi ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarinConverted.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/antimarin/antimarin0311_test.sdf ${SRC_GATHERING_DB_PATH}/antimarin/converting.py
	cd ${PWD}/src && python ${SRC_GATHERING_DB_PATH}/antimarin/converting.py

biofacquim: ${INTERIM_DB_PATH}/biofacquim.tsv.gz
${INTERIM_DB_PATH}/biofacquim.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz ${SRC_GATHERING_DB_PATH}/biofacquim/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/biofacquim/standardizing.R

biofacquim-convert: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/biofacquim/biofacquim.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf ${SRC_GATHERING_DB_PATH}/biofacquim/converting.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/biofacquim/converting.R

biofacquim-download: ${EXTERNAL_DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf
${EXTERNAL_DB_SOURCE_PATH}/biofacquim/BIOFACQUIM_V2.sdf: ${SRC_GATHERING_DB_PATH}/biofacquim/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/biofacquim/downloading.sh

biophytmol: ${INTERIM_DB_PATH}/biophytmol.tsv.gz
${INTERIM_DB_PATH}/biophytmol.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/biophytmol/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/biophytmol/standardizing.R

biophytmol-scrape: ${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/biophytmol/biophytmolScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/biophytmol/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/biophytmol/scraping.R

carotenoiddb: ${INTERIM_DB_PATH}/carotenoiddb.tsv.gz
${INTERIM_DB_PATH}/carotenoiddb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/carotenoiddb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/carotenoiddb/standardizing.R

carotenoiddb-scrape: ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/carotenoiddbScraped.tsv.gz: # ${EXTERNAL_DB_SOURCE_PATH}/carotenoiddb/Carotenoids_InChI_InChIKey.tsv ${SRC_GATHERING_DB_PATH}/carotenoiddb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/carotenoiddb/scraping.R

coconut: ${INTERIM_DB_PATH}/coconut.tsv.gz
${INTERIM_DB_PATH}/coconut.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz ${SRC_GATHERING_DB_PATH}/coconut/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/coconut/standardizing.R

coconut-convert: ${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/coconut/coconutConverted.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/coconut/COCONUT_DB.sdf ${SRC_GATHERING_DB_PATH}/coconut/converting.py
	cd ${PWD}/src && python ${SRC_GATHERING_DB_PATH}/coconut/converting.py

coconut-download: ${EXTERNAL_DB_SOURCE_PATH}/coconut/COCONUT_DB.sdf
${EXTERNAL_DB_SOURCE_PATH}/coconut/COCONUT_DB.sdf: ${SRC_GATHERING_DB_PATH}/coconut/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/coconut/downloading.sh

cyanometdb: ${INTERIM_DB_PATH}/cyanometdb.tsv.gz
${INTERIM_DB_PATH}/cyanometdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/cyanometdb/CyanoMetDB_WR_Feb2021.csv ${SRC_GATHERING_DB_PATH}/cyanometdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/cyanometdb/standardizing.R

cyanometdb-download: ${EXTERNAL_DB_SOURCE_PATH}/cyanometdb/CyanoMetDB_WR_Feb2021.csv
${EXTERNAL_DB_SOURCE_PATH}/cyanometdb/CyanoMetDB_WR_Feb2021.csv: ${SRC_GATHERING_DB_PATH}/cyanometdb/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/cyanometdb/downloading.sh

datawarrior: ${INTERIM_DB_PATH}/datawarrior.tsv.gz
${INTERIM_DB_PATH}/datawarrior.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/datawarrior/NaturalProducts.txt ${SRC_GATHERING_DB_PATH}/datawarrior/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/datawarrior/standardizing.R

datawarrior-download-modified: ${EXTERNAL_DB_SOURCE_PATH}/datawarrior/NaturalProducts.txt
${EXTERNAL_DB_SOURCE_PATH}/datawarrior/NaturalProducts.txt: ${SRC_GATHERING_DB_PATH}/datawarrior/downloading_modified.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/datawarrior/downloading_modified.sh

dianatdb: ${INTERIM_DB_PATH}/dianatdb.tsv.gz
${INTERIM_DB_PATH}/dianatdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/dianatdb/2020_DiaNatDB_336.xlsx ${SRC_GATHERING_DB_PATH}/dianatdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/dianatdb/standardizing.R

dianatdb-download: ${EXTERNAL_DB_SOURCE_PATH}/dianatdb/2020_DiaNatDB_336.xlsx
${EXTERNAL_DB_SOURCE_PATH}/dianatdb/2020_DiaNatDB_336.xlsx: ${SRC_GATHERING_DB_PATH}/dianatdb/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/dianatdb/downloading.sh

dnp: ${INTERIM_DB_PATH}/dnp.tsv.gz
${INTERIM_DB_PATH}/dnp.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/dnp/30_1/full_set.csv ${SRC_GATHERING_DB_PATH}/dnp/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/dnp/standardizing.R

## TODO add dnp extracting step

drduke: ${INTERIM_DB_PATH}/drduke.tsv.gz
${INTERIM_DB_PATH}/drduke.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv ${SRC_GATHERING_DB_PATH}/drduke/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/drduke/standardizing.R

drduke-download: ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv
${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FARMACY_NEW.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/REFERENCES.csv: ${SRC_GATHERING_DB_PATH}/drduke/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/drduke/downloading.sh

foodb: ${INTERIM_DB_PATH}/foodb.tsv.gz
${INTERIM_DB_PATH}/foodb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/CompoundsFlavor_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Compound_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Content.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Flavor.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Reference.csv ${SRC_GATHERING_DB_PATH}/foodb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/foodb/standardizing.R

foodb-download: ## TODO target directory
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/foodb/downloading.sh

foodb-download-modified: ## TODO target multiple files
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/foodb/downloading_modified.sh

inflamnat: ${INTERIM_DB_PATH}/inflamnat.tsv.gz
${INTERIM_DB_PATH}/inflamnat.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx ${SRC_GATHERING_DB_PATH}/inflamnat/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/inflamnat/standardizing.R

inflamnat-download: ${EXTERNAL_DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx
${EXTERNAL_DB_SOURCE_PATH}/inflamnat/ci8b00560_si_001.xlsx: ${SRC_GATHERING_DB_PATH}/inflamnat/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/inflamnat/downloading.sh

knapsack: ${INTERIM_DB_PATH}/knapsack.tsv.gz
${INTERIM_DB_PATH}/knapsack.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/knapsack/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/knapsack/standardizing.R

knapsack-scrape: ${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/knapsack/knapsackScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/knapsack/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/knapsack/scraping.R

metabolights: ${INTERIM_DB_PATH}/metabolights.tsv.gz
${INTERIM_DB_PATH}/metabolights.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml ${SRC_GATHERING_DB_PATH}/metabolights/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/metabolights/standardizing.R

metabolights-download: ${EXTERNAL_DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml
${EXTERNAL_DB_SOURCE_PATH}/metabolights/eb-eye_metabolights_complete.xml: ${SRC_GATHERING_DB_PATH}/metabolights/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/metabolights/downloading.sh

mibig: ${INTERIM_DB_PATH}/mibig.tsv.gz
${INTERIM_DB_PATH}/mibig.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/mibig/mibig_json_2.0.zip ${SRC_GATHERING_DB_PATH}/mibig/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/mibig/standardizing.R

mibig-download: ## TODO target directory
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/mibig/downloading.sh

mitishamba: ${INTERIM_DB_PATH}/mitishamba.tsv.gz
${INTERIM_DB_PATH}/mitishamba.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/mitishamba/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/mitishamba/standardizing.R

mitishamba-scrape: ${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/mitishamba/mitishambaScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/mitishamba/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/mitishamba/scraping.R

nanpdb: ${INTERIM_DB_PATH}/nanpdb.tsv.gz
${INTERIM_DB_PATH}/nanpdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/nanpdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/nanpdb/standardizing.R

nanpdb-scrape: ${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/nanpdb/nanpdbScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/nanpdb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/nanpdb/scraping.R

napralert: ${INTERIM_DB_PATH}/napralert.tsv.gz
${INTERIM_DB_PATH}/napralert.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert_matched_final_unified.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert_matched_final_unified.tsv.gz ${EXTERNAL_DB_SOURCE_PATH}/napralert/napralert.tsv.gz ${SRC_GATHERING_DB_PATH}/napralert/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/napralert/standardizing.R

napralert-download: ## TODO target multiple files
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/napralert/downloading.sh

npass: ${INTERIM_DB_PATH}/npass.tsv.gz
${INTERIM_DB_PATH}/npass.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_properties.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt ${EXTERNAL_DB_SOURCE_PATH}/npass/NPASSv1.0_download_naturalProducts_species_pair.txt ${SRC_GATHERING_DB_PATH}/npass/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/npass/standardizing.R

npass-download: ## TODO target multiple files
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/npass/downloading.sh

npatlas: ${INTERIM_DB_PATH}/npatlas.tsv.gz
${INTERIM_DB_PATH}/npatlas.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npatlas/NPAtlas_download.tsv ${SRC_GATHERING_DB_PATH}/npatlas/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/npatlas/standardizing.R

npatlas-download: ${EXTERNAL_DB_SOURCE_PATH}/npatlas/NPAtlas_download.tsv
${EXTERNAL_DB_SOURCE_PATH}/npatlas/NPAtlas_download.tsv: ${SRC_GATHERING_DB_PATH}/npatlas/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/npatlas/downloading.sh

npcare: ${INTERIM_DB_PATH}/npcare.tsv.gz
${INTERIM_DB_PATH}/npcare.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npcare/npcare.zip ${SRC_GATHERING_DB_PATH}/npcare/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/npcare/standardizing.R

npcare-download: ${EXTERNAL_DB_SOURCE_PATH}/npcare/npcare.zip
${EXTERNAL_DB_SOURCE_PATH}/npcare/npcare.zip: ${SRC_GATHERING_DB_PATH}/npcare/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/npcare/downloading.sh

npedia: ${INTERIM_DB_PATH}/npedia.tsv.gz
${INTERIM_DB_PATH}/npedia.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/npedia/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/npedia/standardizing.R

npedia-scrape: ${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/npedia/npediaScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/npedia/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/npedia/scraping.R

nubbe: ${INTERIM_DB_PATH}/nubbe.tsv.gz
${INTERIM_DB_PATH}/nubbe.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/nubbe/nubbe_raw/*.xml) ${SRC_GATHERING_DB_PATH}/nubbe/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/nubbe/standardizing.R

# nubbe-scrape TODO

pamdb: ${INTERIM_DB_PATH}/pamdb.tsv.gz
${INTERIM_DB_PATH}/pamdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/pamdb/PaMet.xlsx ${SRC_GATHERING_DB_PATH}/pamdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/pamdb/standardizing.R

pamdb-download: ${EXTERNAL_DB_SOURCE_PATH}/pamdb/PaMet.xlsx
${EXTERNAL_DB_SOURCE_PATH}/pamdb/PaMet.xlsx: ${SRC_GATHERING_DB_PATH}/pamdb/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/pamdb/downloading.sh

phenolexplorer: ${INTERIM_DB_PATH}/phenolexplorer.tsv.gz
${INTERIM_DB_PATH}/phenolexplorer.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/compounds-structures.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/compounds.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/publications.csv ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/composition-data.xlsx ${SRC_GATHERING_DB_PATH}/phenolexplorer/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/phenolexplorer/standardizing.R

phenolexplorer-download: ## TODO target multiple files
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/phenolexplorer/downloading.sh

phytohub: ${INTERIM_DB_PATH}/phytohub.tsv.gz
${INTERIM_DB_PATH}/phytohub.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/phytohub/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/phytohub/standardizing.R

phytohub-scrape: ${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/phytohub/phytohubScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/phytohub/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/phytohub/scraping.R

procardb: ${INTERIM_DB_PATH}/procardb.tsv.gz
${INTERIM_DB_PATH}/procardb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/procardb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/procardb/standardizing.R

procardb-scrape: ${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/procardb/procardbScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/procardb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/procardb/scraping.R

respect: ${INTERIM_DB_PATH}/respect.tsv.gz
${INTERIM_DB_PATH}/respect.tsv.gz: $(wildcard ${EXTERNAL_DB_SOURCE_PATH}/respect/*.txt) ${SRC_GATHERING_DB_PATH}/respect/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/respect/standardizing.R

respect-download: ## TODO target directory
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/respect/downloading.sh

sancdb: ${INTERIM_DB_PATH}/sancdb.tsv.gz
${INTERIM_DB_PATH}/sancdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/sancdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/sancdb/standardizing.R

sancdb-scrape: ${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/sancdb/sancdbScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/sancdb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/sancdb/scraping.R

streptomedb: ${INTERIM_DB_PATH}/streptomedb.tsv.gz
${INTERIM_DB_PATH}/streptomedb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz ${SRC_GATHERING_DB_PATH}/streptomedb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/streptomedb/standardizing.R

streptomedb-convert: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.sdf ${SRC_GATHERING_DB_PATH}/streptomedb/converting.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/streptomedb/converting.R

streptomedb-download: ${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.sdf
${EXTERNAL_DB_SOURCE_PATH}/streptomedb/streptomedb.sdf: ${SRC_GATHERING_DB_PATH}/streptomedb/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/streptomedb/downloading.sh

swmd: ${INTERIM_DB_PATH}/swmd.tsv.gz
${INTERIM_DB_PATH}/swmd.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/swmd/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/swmd/standardizing.R

swmd-scrape: ${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/swmd/swmdScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/swmd/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/swmd/scraping.R

tmdb: ${INTERIM_DB_PATH}/tmdb.tsv.gz
${INTERIM_DB_PATH}/tmdb.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/tmdb/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/tmdb/standardizing.R

tmdb-scrape: ${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/tmdb/tmdbScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/tmdb/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/tmdb/scraping.R

tmmc: ${INTERIM_DB_PATH}/tmmc.tsv.gz
${INTERIM_DB_PATH}/tmmc.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx ${SRC_GATHERING_DB_PATH}/tmmc/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/tmmc/standardizing.R

tmmc-download: ${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx
${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx: ${SRC_GATHERING_DB_PATH}/tmmc/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/tmmc/downloading.sh

tppt: ${INTERIM_DB_PATH}/tppt.tsv.gz
${INTERIM_DB_PATH}/tppt.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/tppt/TPPT_database.xlsx ${SRC_GATHERING_DB_PATH}/tppt/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/tppt/standardizing.R

tppt-download: ${EXTERNAL_DB_SOURCE_PATH}/tppt/TPPT_database.xlsx
${EXTERNAL_DB_SOURCE_PATH}/tppt/TPPT_database.xlsx: ${SRC_GATHERING_DB_PATH}/tppt/downloading.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/tppt/downloading.sh

unpd: ${INTERIM_DB_PATH}/unpd.tsv.gz
${INTERIM_DB_PATH}/unpd.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz ${SRC_GATHERING_DB_PATH}/unpd/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/unpd/standardizing.R

unpd-download: ## TODO target multiple files
	cd ${PWD}/src && bash ${SRC_GATHERING_DB_PATH}/unpd/downloading.sh

unpd-integrate: ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/unpd/unpdIntegrated.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/unpd/UNPD_DB.csv.zip ${EXTERNAL_DB_SOURCE_PATH}/unpd/unpd_final.csv.zip ${SRC_GATHERING_DB_PATH}/unpd/integrating.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/unpd/integrating.R

wakankensaku: ${INTERIM_DB_PATH}/wakankensaku.tsv.gz
${INTERIM_DB_PATH}/wakankensaku.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz ${SRC_GATHERING_DB_PATH}/wakankensaku/standardizing.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/wakankensaku/standardizing.R

wakankensaku-scrape: ${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz
${EXTERNAL_DB_SOURCE_PATH}/wakankensaku/wakankensakuScraped.tsv.gz: # ${SRC_GATHERING_DB_PATH}/wakankensaku/scraping.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_DB_PATH}/wakankensaku/scraping.R
