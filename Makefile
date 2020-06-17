.PHONY: help docker-build docker-bash databases afrotryp

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

databases: afrotryp

afrotryp: ${DATA_PATH}/interim/db/afrotryp.tsv.zip

${DATA_PATH}/interim/db/afrotryp.tsv.zip: ${DATA_PATH}/external/dbSource/afrotryp/afrotryp.tsv.zip
	cd src &&	Rscript 1_gathering/db/afrotryp/standardizing.R
