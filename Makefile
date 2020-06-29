include paths.mk

.PHONY: help docker-build docker-bash databases
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
	docker build -t onpdb-environment --build-arg USER_ID=$(shell id -u) --build-arg GROUP_ID=$(shell id -g) .

docker-bash:
	docker run -it --rm -v $$PWD:/srv/onpdb onpdb-environment 

databases:
	make -C ${SRC_GATHERING_PATH} databases

databases-reconvert: ${DATABASES_RECONVERT}
	make -C ${SRC_GATHERING_PATH} databases-reconvert

databases-reintegrate: ${DATABASES_REINTEGRATE}
	make -C ${SRC_GATHERING_PATH} databases-reintegrate

databases-rescrape: ${DATABASES_RESCRAPE}
	make -C ${SRC_GATHERING_PATH} databases-rescrape

curating: curating-integrating curating-editing curating-editing-chemo

curating-integrating: ${DATABASES}
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

curating-editing-chemo-smiles: ${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureSmiles_min.tsv.zip
${INTERIM_TABLE_TRANSLATED_PATH}/translatedStructureSmiles_min.tsv.zip: ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/smiles_min.py
${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH}/smiles_min.py: ${INTERIM_TABLE_ORIGINAL_PATH}/originalStructureSmiles.tsv.zip
	cd src && python 2_curating/2_editing/chemo/subscripts/1_translating/smiles_min.py
