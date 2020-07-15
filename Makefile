include config.mk
include paths.mk

.PHONY: help docker-build docker-bash 
.PHONY: gathering-full gathering-databases-full gathering-databases gathering-databases-reconvert gathering-databases-reintegrate gathering-databases-rescrape gathering-translation-full gathering-translation-common gathering-translation-tcm
.PHONY: curating curating-1-integrating curating-editing curating-3-integrating
.PHONY: curating-editing-organism curating-editing-organism-cleaning-original curating-editing-organism-translating curating-editing-organism-cleaning-translated curating-editing-organism-cleaning-taxonomy
.PHONY: curating-editing-reference curating-editing-reference-translating curating-editing-reference-translating-doi curating-editing-reference-translating-pubmed curating-editing-reference-translating-title curating-editing-reference-translating-unsplit curating-editing-reference-integrating curating-editing-reference-cleaning
.PHONY: curating-editing-structure curating-editing-structure-translating curating-editing-structure-translating-name curating-editing-structure-translating-smiles curating-editing-structure-integrating curating-editing-structure-sanitizing
.PRECIOUS: %.tsv %.zip %.json %.gz

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

gathering-full: gathering-databases-full gathering-translation-full

gathering-databases-full: gathering-databases gathering-databases-reconvert gathering-databases-reintegrate gathering-databases-rescrape

gathering-databases:
	make	-C	${SRC_GATHERING_DB_PATH}	gathering-databases

gathering-databases-reconvert: ${DATABASES_RECONVERT}
	make	-C	${SRC_GATHERING_DB_PATH}	gathering-databases-reconvert

gathering-databases-reintegrate: ${DATABASES_REINTEGRATE}
	make	-C	${SRC_GATHERING_DB_PATH}	gathering-databases-reintegrate

gathering-databases-rescrape: ${DATABASES_RESCRAPE}
	make	-C	${SRC_GATHERING_DB_PATH}	gathering-databases-rescrape

gathering-translation-full: gathering-translation-common gathering-translation-tcm

gathering-translation-common:
	make	-C	${SRC_GATHERING_TRANSLATION_PATH}	gathering-translation-common

gathering-translation-tcm:
	make	-C	${SRC_GATHERING_TRANSLATION_PATH}	gathering-translation-tcm

curating: curating-1-integrating curating-editing curating-3-integrating

curating-1-integrating: ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz: ${DATABASES} paths.mk ${SRC_PATH}/paths.R ${SRC_CURATING_1_INTEGRATING_PATH}/integrating.R
	cd	src	&&	Rscript	${SRC_CURATING_1_INTEGRATING_PATH}/integrating.R

curating-editing: curating-editing-organism curating-editing-reference curating-editing-structure

curating-editing-organism: curating-editing-organism-cleaning-original curating-editing-organism-translating curating-editing-organism-cleaning-translated curating-editing-organism-cleaning-taxonomy

curating-editing-organism-cleaning-original: ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/original.tsv.gz
${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/original.tsv.gz: $(wildcard ${INTERIM_TABLE_ORIGINAL_ORGANISM_PATH}/*.tsv) ${INTERIM_DICTIONARY_PATH}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/1_cleaningOriginal.R
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/1_cleaningOriginal.R

curating-editing-organism-translating: ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/interim.tsv.gz
${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH_KT}/build/libs/shadow.jar: ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH_KT}/build.gradle.kts $(wildcard ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH_KT}/src/main/kotlin/*.kt)
	./gradlew castShadows
${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/interim.tsv.gz: ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_DICTIONARY_PATH}/common/black.tsv ${INTERIM_DICTIONARY_PATH}/common/manualSubtraction.tsv ${INTERIM_DICTIONARY_PATH}/common/names.tsv.gz ${INTERIM_DICTIONARY_PATH}/tcm/names.tsv.gz ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/2_translating_organism_kotlin/build/libs/shadow.jar
	@java -jar ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH_KT}/build/libs/shadow.jar ${DATA_PATH} ${FULL}

curating-editing-organism-cleaning-translated: curating-editing-organism-translating ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/translated.tsv.gz
${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/translated.tsv.gz: $(wildcard ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/translated/*.json) ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_DICTIONARY_PATH}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/3_cleaningTranslated.R
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/3_cleaningTranslated.R

curating-editing-organism-cleaning-taxonomy: ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/cleaned.tsv.gz
${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/cleaned.tsv.gz: ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/translated.tsv.gz ${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/4_cleaningTaxonomy.R
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_ORGANISM_SUBSCRIPTS_PATH}/4_cleaningTaxonomy.R

curating-editing-reference: curating-editing-reference-translating curating-editing-reference-integrating curating-editing-reference-cleaning

curating-editing-reference-translating: curating-editing-reference-translating-doi curating-editing-reference-translating-pubmed curating-editing-reference-translating-title curating-editing-reference-translating-unsplit

curating-editing-reference-translating-doi: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/doi.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/doi.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/doi.R

curating-editing-reference-translating-pubmed: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/pubmed.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/pubmed.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/pubmed.R

curating-editing-reference-translating-title: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/title.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/title.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/title.R

curating-editing-reference-translating-unsplit: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/unsplit.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/unsplit.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/unsplit.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/unsplit.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_TRANSLATING_PATH}/unsplit.R

curating-editing-reference-integrating: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_PATH}/2_integrating.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz  ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/unsplit.tsv.gz ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/full.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_PATH}/2_integrating.R

curating-editing-reference-cleaning: ${INTERIM_TABLE_CLEANED_REFERENCE_PATH}/cleaned.tsv.gz
${INTERIM_TABLE_CLEANED_REFERENCE_PATH}/cleaned.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_PATH}/3_cleaning.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_REFERENCE_SUBSCRIPTS_PATH}/3_cleaning.R

curating-editing-structure: curating-editing-structure-translating curating-editing-structure-integrating curating-editing-structure-sanitizing

curating-editing-structure-translating: curating-editing-structure-translating-name curating-editing-structure-translating-smiles

curating-editing-structure-translating-name: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_TRANSLATING_PATH}/names.R ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/nominal.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_TRANSLATING_PATH}/names.R

curating-editing-structure-translating-smiles: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_TRANSLATING_PATH}/smiles.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/smiles.tsv.gz
	cd	src	&&	python	${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_TRANSLATING_PATH}/smiles.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz structureOriginalSmiles

curating-editing-structure-integrating: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_INTEGRATING_PATH}/integrating.R ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_INTEGRATING_PATH}/integrating.R

curating-editing-structure-sanitizing: ${INTERIM_TABLE_CLEANED_STRUCTURE_PATH}/cleaned.tsv.gz
${INTERIM_TABLE_CLEANED_STRUCTURE_PATH}/cleaned.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_CLEANINGANDENRICHING_PATH}/chemosanitizer.py ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz
	cd	src	&&	python	${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_CLEANINGANDENRICHING_PATH}/chemosanitizer.py ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz ${INTERIM_TABLE_CLEANED_STRUCTURE_PATH}/cleaned.tsv.gz structureTranslated 8

curating-3-integrating:	${INTERIM_TABLE_CLEANED_PATH}/table.tsv.gz
${INTERIM_TABLE_CLEANED_PATH}/table.tsv.gz: ${SRC_CURATING_3_INTEGRATING_PATH}/integrating.R ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_CLEANED_STRUCTURE_PATH}/cleaned.tsv.gz ${INTERIM_TABLE_CLEANED_ORGANISM_PATH}/cleaned.tsv.gz ${INTERIM_TABLE_CLEANED_REFERENCE_PATH}/cleaned.tsv.gz
	cd	src	&&	Rscript	${SRC_CURATING_3_INTEGRATING_PATH}/integrating.R
