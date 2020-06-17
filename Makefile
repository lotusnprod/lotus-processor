DATA_PATH ?= ${PWD}/data

INTERIM_PATH = ${DATA_PATH}/interim/db
SOURCE_PATH = ${DATA_PATH}/external/dbSource

AFROTRYP_SOURCE_PATH = ${SOURCE_PATH}/afrotryp
ALKAMID_SOURCE_PATH = ${SOURCE_PATH}/alkamid


.PHONY: help docker-build docker-bash databases afrotryp alkamid alkamid-rescrape

help:
	@echo "Builder"
	@echo "-------"
	@echo ""
	@echo "docker-build: build the docker image (with no data)"
	@echo "docker-bash: run a shell into the docker image"

docker-build:
	docker build -t onpdb-environment .

docker-bash:
	docker run -it --rm -v $$PWD:/srv/onpdb onpdb-environment bash

databases: afrotryp alkamid

databases-rescrape: alkamid-rescrape

afrotryp: ${INTERIM_PATH}/afrotryp.tsv.zip

${DATA_PATH}/interim/db/afrotryp.tsv.zip: ${DATA_PATH}/external/dbSource/afrotryp/afrotryp.tsv.zip
	cd src &&	Rscript 1_gathering/db/afrotryp/standardizing.R

alkamid: ${INTERIM_PATH}/alkamid.tsv.zip

${DATA_PATH}/interim/db/alkamid.tsv.zip: ${ALKAMID_SOURCE_PATH}/alkamidRefScraped.tsv.zip ${ALKAMID_SOURCE_PATH}/alkamidRefScraped.tsv.zip
	cd src &&	Rscript 1_gathering/db/alkamid/standardizing.R

alkamid-rescrape:
	cd src && Rscript 1_gathering/db/alkamid/scraping.R
