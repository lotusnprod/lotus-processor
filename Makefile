include config.mk
include paths.mk
include src/1_gathering/db/Makefile
include src/1_gathering/translation/Makefile

.PHONY: help docker-build docker-bash tests
.PHONY: gathering-full gathering-databases-full gathering-databases-full-quick gathering-databases gathering-databases-convert gathering-databases-download gathering-databases-integrate gathering-databases-scrape gathering-translation-full gathering-translation-common gathering-translation-tcm
.PHONY: curating curating-1-integrating curating-editing curating-3-integrating
.PHONY: curating-editing-structure curating-editing-structure-translating curating-editing-structure-translating-name curating-editing-structure-translating-smiles curating-editing-structure-integrating curating-editing-structure-sanitizing curating-editing-structure-naming curating-editing-structure-classifying
.PHONY: curating-editing-organism curating-editing-organism-processing-original curating-editing-organism-translating curating-editing-organism-processing-translated curating-editing-organism-processing-taxonomy
.PHONY: curating-editing-reference curating-editing-reference-translating curating-editing-reference-translating-doi curating-editing-reference-translating-pubmed curating-editing-reference-translating-title curating-editing-reference-translating-split curating-editing-reference-translating-publishingDetails curating-editing-reference-translating-original curating-editing-reference-integrating curating-editing-reference-processing
.PHONY: curating-and-analysing analysing analysing-sampling analysing-validating analysing-metrics analysing-examples
.PHONY: processing-organism-interim
.PHONY: curating-and-analysing-and-visalizing visualizing visualizing-alluvial visualizing-chord visualizing-tree visualizing-upset visualizing-distribution
.PHONY: get-gnfinder get-gnverifier get-opsin get-bins
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
	@echo "get-bins: download gnfinder, gnverify and opsin"
	@echo ""
	@echo "gathering-full: Run the 1_gathering scripts"
	@echo "curating: Run the 2_curating scripts"
	@echo "analysing: Run the 3_analysing scripts"
	@echo "visualizing: Run the 4_visualizing scripts"

docker-build:
	docker build -t onpdb-environment --build-arg USER_ID=$(shell id -u) --build-arg GROUP_ID=$(shell id -g) .

docker-bash:
	docker run -it --rm -v $$PWD:/srv/onpdb onpdb-environment

get-bins: get-gnfinder get-gnverifier get-opsin

get-gnfinder: bin/gnfinder
get-gnverifier: bin/gnverifier
get-opsin: bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar

bin/gnfinder: ${BIN_PATH}/gnfinder
${BIN_PATH}/gnfinder: config.mk
	mkdir -p bin
	curl -L https://github.com/gnames/gnfinder/releases/download/${GNFINDER_VERSION}/gnfinder-${GNFINDER_VERSION}-${PLATFORM}.tar.gz | tar xOz gnfinder > bin/gnfinder
	chmod +x bin/gnfinder

bin/gnverifier: ${BIN_PATH}/gnverifier
${BIN_PATH}/gnverifier: config.mk
	mkdir -p bin
	curl -L https://github.com/gnames/gnverifier/releases/download/${GNVERIFIER_VERSION}/gnverifier-${GNVERIFIER_VERSION}-${PLATFORM}.tar.gz | tar xOz gnverifier > bin/gnverifier
	chmod +x bin/gnverifier

bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar: ${BIN_PATH}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar
${BIN_PATH}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar: config.mk
	mkdir -p bin
	curl -L https://github.com/dan2097/opsin/releases/download/${OPSIN_VERSION}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar > bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar
	chmod +x bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar

tests:
	cd src && Rscript ${TESTS_PATH}/tests.R 

gathering-full: gathering-custom-dictionaries gathering-databases-full gathering-translation-full

gathering-full-quick: gathering-custom-dictionaries gathering-databases-full-quick gathering-translation-full

gathering-custom-dictionaries: 
	cd src && bash ${SRC_GATHERING_PATH}/gathering_custom_dictionaries.sh

gathering-databases-full: gathering-databases-scrape gathering-databases-download-modified gathering-databases-convert gathering-databases-integrate gathering-databases

gathering-databases-full-quick: gathering-databases-download gathering-databases-download-modified gathering-databases-convert gathering-databases-integrate gathering-databases

gathering-databases:
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases

gathering-databases-convert: ${DATABASES_CONVERT}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-convert

gathering-databases-download: ${DATABASES_DOWNLOAD}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-download

gathering-databases-download-modified: ${DATABASES_DOWNLOAD}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-download-modified

gathering-databases-integrate: ${DATABASES_INTEGRATE}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-integrate

gathering-databases-scrape: ${DATABASES_SCRAPE}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-scrape

gathering-translation-full: gathering-translation-common gathering-translation-tcm

gathering-translation-common:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-common

gathering-translation-tcm:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-tcm

curating-and-analysing-and-visualizing: curating analysing visualizing

curating-and-analysing: curating analysing

curating: get-bins curating-1-integrating curating-editing curating-3-integrating

curating-1-integrating: ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz: ${DATABASES} paths.mk ${SRC_PATH}/paths.R ${SRC_CURATING_PATH}/1_integrating.R
	cd src && Rscript ${SRC_CURATING_PATH}/1_integrating.R

curating-editing: curating-editing-structure curating-editing-organism curating-editing-reference

curating-editing-structure: curating-editing-structure-translating curating-editing-structure-integrating curating-editing-structure-sanitizing curating-editing-structure-stereocounting curating-editing-structure-naming curating-editing-structure-classifying

curating-editing-structure-translating: curating-editing-structure-translating-name curating-editing-structure-translating-smiles

curating-editing-structure-translating-name: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/names.R ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/nominal.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/names.R

curating-editing-structure-translating-smiles: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/smiles.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/smiles.tsv.gz
	cd src && python ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/smiles.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz structureOriginal_smiles

curating-editing-structure-integrating: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_PATH}/2_integrating.R ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_STRUCTURE_PATH}/2_integrating.R

curating-editing-structure-sanitizing: ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz
${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_PROCESSINGANDENRICHING_PATH}/chemosanitizer.py ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz
	cd src && python ${SRC_CURATING_EDITING_STRUCTURE_PROCESSINGANDENRICHING_PATH}/chemosanitizer.py ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz structureTranslated 8

curating-editing-structure-stereocounting: ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/counted.tsv.gz
${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/counted.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_PROCESSINGANDENRICHING_PATH}/stereocounter.py ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz
	cd src && python ${SRC_CURATING_EDITING_STRUCTURE_PROCESSINGANDENRICHING_PATH}/stereocounter.py ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/counted.tsv.gz smilesSanitized

curating-editing-structure-naming: # add dependent files
	cd src && Rscript ${SRC_CURATING_EDITING_STRUCTURE_ENRICHING_PATH}/naming.R

curating-editing-structure-classifying: # add dependent files
	cd src && Rscript ${SRC_CURATING_EDITING_STRUCTURE_ENRICHING_PATH}/np_classifier.R

curating-editing-organism: curating-editing-organism-processing-original processing-organism-interim curating-editing-organism-translating curating-editing-organism-processing-translated curating-editing-organism-processing-taxonomy

curating-editing-organism-processing-original: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz: $(wildcard ${INTERIM_TABLE_ORIGINAL_ORGANISM_PATH}/^[0-9]{6}.tsv) ${INTERIM_DICTIONARY_PATH_FIX}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_PATH}/1_processingOriginal.R
	cd src && Rscript ${SRC_CURATING_EDITING_ORGANISM_PATH}/1_processingOriginal.R

processing-organism-interim:
	-rm edit ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz

curating-editing-organism-translating: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz
${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build/libs/shadow.jar: ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build.gradle.kts $(wildcard ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/src/main/kotlin/*.kt)
	./gradlew castShadows
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_DICTIONARY_PATH_FIX}/common/deny.tsv ${INTERIM_DICTIONARY_PATH_FIX}/common/manualSubtraction.tsv ${INTERIM_DICTIONARY_PATH_FIX}/common/names.tsv.gz ${INTERIM_DICTIONARY_PATH_FIX}/tcm/names.tsv.gz ${SRC_CURATING_EDITING_ORGANISM_PATH}/2_translating_organism_kotlin/build/libs/shadow.jar
	@java -jar ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build/libs/shadow.jar ${DATA_PATH} ${MODE}

curating-editing-organism-processing-translated: curating-editing-organism-translating ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated.tsv.gz
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated.tsv.gz: $(wildcard ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated/*.json) ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_DICTIONARY_PATH_FIX}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_PATH}/3_processingTranslated.R
	cd src && Rscript ${SRC_CURATING_EDITING_ORGANISM_PATH}/3_processingTranslated.R

curating-editing-organism-processing-taxonomy: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/processed.tsv.gz
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/processed.tsv.gz: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated.tsv.gz ${SRC_CURATING_EDITING_ORGANISM_PATH}/4_processingTaxonomy.R
	cd src && Rscript ${SRC_CURATING_EDITING_ORGANISM_PATH}/4_processingTaxonomy.R

curating-editing-reference: curating-editing-reference-translating curating-editing-reference-integrating curating-editing-reference-processing

curating-editing-reference-translating: curating-editing-reference-translating-doi curating-editing-reference-translating-pubmed curating-editing-reference-translating-title curating-editing-reference-translating-split curating-editing-reference-translating-publishingDetails curating-editing-reference-translating-original

curating-editing-reference-translating-doi: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/doi.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/doi.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/doi.R

curating-editing-reference-translating-pubmed: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/pubmed.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/pubmed.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/pubmed.R

curating-editing-reference-translating-title: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/title.R $(wildcard ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/title/*.tsv.gz)
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/title.R

curating-editing-reference-translating-split: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/split.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/split.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/split.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/split.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/split.R

curating-editing-reference-translating-publishingDetails: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/publishingDetails.R ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/publishingDetails.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/publishingDetails.R

curating-editing-reference-translating-original: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/original.R $(wildcard ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/original/*.tsv.gz)
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/original.R

curating-editing-reference-integrating: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_PATH}/2_integrating.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/split.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/full.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_PATH}/2_integrating.R

curating-editing-reference-processing: ${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz
${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_PATH}/3_processing.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz ${EXTERNAL_TRANSLATION_SOURCE_PATH}/pubmed/PMC-ids.csv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_PATH}/3_processing.R

curating-3-integrating: ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz: ${SRC_CURATING_PATH}/3_integrating.R ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/smiles.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/processed.tsv.gz ${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz
	cd src && Rscript ${SRC_CURATING_PATH}/3_integrating.R

analysing: get-bins analysing-sampling analysing-validating analysing-metrics analysing-examples

analysing-sampling: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/1_sampling.R

analysing-validating: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/2_validating.R

analysing-metrics: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/3_metrics.R

analysing-examples: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz ${SRC_ANALYSING_PATH}/examples.R
	cd src && Rscript ${SRC_ANALYSING_PATH}/examples.R

visualizing: visualizing-alluvial visualizing-chord visualizing-tree visualizing-upset visualizing-distribution

visualizing-alluvial:
	cd src && Rscript ${SRC_VISUALIZING_PATH}/plot_alluvial.R

visualizing-chord:
	cd src && Rscript ${SRC_VISUALIZING_PATH}/plot_chordDiagrams.R

visualizing-tree:
	cd src && Rscript ${SRC_VISUALIZING_PATH}/plot_magicTree.R

visualizing-upset:
	cd src && Rscript ${SRC_VISUALIZING_PATH}/plot_upset.R

visualizing-distribution:
	cd src && Rscript ${SRC_VISUALIZING_PATH}/plot_distribution.R
