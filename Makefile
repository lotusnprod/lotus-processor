include config.mk
include paths.mk
include ${SRC_GATHERING_DB_PATH}/Makefile
include ${SRC_GATHERING_TAXONOMY_PATH}/Makefile
include ${SRC_GATHERING_TRANSLATION_PATH}/Makefile

.PHONY: help docker-build docker-bash tests
.PHONY: gathering-full gathering-quick gathering-full-hard 
.PHONY: gathering-databases-full gathering-databases-full-hard gathering-databases-quick 
.PHONY: gathering-databases gathering-databases-download gathering-databases-download-modified 
.PHONY: gathering-databases-scrape gathering-databases-accessible gathering-databases-semi gathering-databases-closed
.PHONY: gathering-custom-dictionaries gathering-pmcid gathering-gbif gathering-chinese-board gathering-validation
.PHONY: gathering-translation-tcmid gathering-translation-full gathering-translation-common gathering-translation-common-quick gathering-translation-quick gathering-translation-tcm gathering-translation-tcm-quick
.PHONY: gathering-taxonomy-otl gathering-taxonomy-npclassifier gathering-taxonomy-classyfire gathering-taxonomy-full
.PHONY: curating curating-1-integrating curating-editing curating-3-integrating
.PHONY: curating-editing-structure curating-editing-structure-translating curating-editing-structure-translating-name curating-editing-structure-translating-inchi curating-editing-structure-integrating curating-editing-structure-sanitizing curating-editing-structure-naming curating-editing-structure-classifying
.PHONY: curating-editing-organism curating-editing-organism-processing-original curating-editing-organism-translating curating-editing-organism-processing-translated curating-editing-organism-processing-taxonomy
.PHONY: curating-editing-reference curating-editing-reference-translating curating-editing-reference-translating-doi curating-editing-reference-translating-pubmed curating-editing-reference-translating-title curating-editing-reference-translating-split curating-editing-reference-translating-publishingDetails curating-editing-reference-translating-original curating-editing-reference-integrating curating-editing-reference-processing
.PHONY: curating-and-analyzing analyzing analyzing-sampling analyzing-validating analyzing-metrics analyzing-examples
.PHONY: processing-organism-interim
.PHONY: curating-and-analyzing-and-visalizing visualizing visualizing-alluvial visualizing-chord visualizing-tree visualizing-upset visualizing-distribution
.PHONY: get-gnfinder get-gnverifier get-opsin get-bins
.PHONY: manual-entry lotus lotus-quick
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
	@echo "analyzing: Run the 3_analyzing scripts"
	@echo "visualizing: Run the 4_visualizing scripts"

docker-build:
	docker build -t onpdb-environment --build-arg USER_ID=$(shell id -u) --build-arg GROUP_ID=$(shell id -g) .

docker-bash:
	docker run -it --rm -v $$PWD:/srv/onpdb onpdb-environment

get-bins: get-gnfinder get-gnverifier get-opsin

get-gnfinder: ${BIN_PATH}/gnfinder
get-gnverifier: ${BIN_PATH}/gnverifier
get-opsin: ${BIN_PATH}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar

bin/gnfinder: ${BIN_PATH}/gnfinder
${BIN_PATH}/gnfinder: config.mk
	mkdir -p bin
ifeq (${PLATFORM},$(filter ${PLATFORM},mac linux))
	curl -L https://github.com/gnames/gnfinder/releases/download/${GNFINDER_VERSION}/gnfinder-${GNFINDER_VERSION}-${PLATFORM}.tar.gz | tar xOz gnfinder > bin/gnfinder
	chmod +x bin/gnfinder
else
	curl -L https://github.com/gnames/gnfinder/releases/download/${GNFINDER_VERSION}/gnfinder-${GNFINDER_VERSION}-win-64.zip > bin/gnfinder.zip
	unzip -o bin/gnfinder.zip -d bin/
	rm bin/gnfinder.zip
endif

bin/gnverifier: ${BIN_PATH}/gnverifier
${BIN_PATH}/gnverifier: config.mk
	mkdir -p bin
ifeq (${PLATFORM},$(filter ${PLATFORM},mac linux))
	curl -L https://github.com/gnames/gnverifier/releases/download/${GNVERIFIER_VERSION}/gnverifier-${GNVERIFIER_VERSION}-${PLATFORM}.tar.gz | tar xOz gnverifier > bin/gnverifier
	chmod +x bin/gnverifier
else
	curl -L https://github.com/gnames/gnverifier/releases/download/${GNVERIFIER_VERSION}/gnverifier-${GNVERIFIER_VERSION}-win-64.zip > bin/gnverifier.zip
	unzip -o bin/gnverifier.zip -d bin/
	rm bin/gnverifier.zip
endif

bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar: ${BIN_PATH}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar
${BIN_PATH}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar: config.mk
	mkdir -p bin
	curl -L https://github.com/dan2097/opsin/releases/download/${OPSIN_VERSION}/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar > bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar
	chmod +x bin/opsin-${OPSIN_VERSION}-jar-with-dependencies.jar

tests:
	cd src && Rscript ${TESTS_PATH}/tests.R 

gathering-quick: gathering-custom-dictionaries gathering-translation-quick gathering-taxonomy-quick

gathering-full: gathering-custom-dictionaries gathering-databases-full gathering-translation-full gathering-taxonomy-full

gathering-full-hard: gathering-custom-dictionaries gathering-databases-full-hard gathering-translation-full gathering-taxonomy-full

gathering-custom-dictionaries: 
	cd src && bash ${SRC_GATHERING_PATH}/dictionary/gathering_custom_dictionaries.sh

gathering-validation: 
	cd src && bash ${SRC_GATHERING_VALIDATION_PATH}/gathering_validation.sh

gathering-databases-quick: gathering-databases-download gathering-databases-download-modified gathering-databases-accessible

gathering-databases-full: gathering-databases-scrape gathering-databases-quick

gathering-databases-full-hard: gathering-databases-semi gathering-databases-closed gathering-databases-full

gathering-databases-accessible:
	mkdir -p ${INTERIM_DB_PATH}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-accessible

gathering-databases-semi:
	mkdir -p ${INTERIM_DB_PATH}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-semi

gathering-databases-closed:
	mkdir -p ${INTERIM_DB_PATH}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-closed

gathering-databases-download: ${DATABASES_DOWNLOAD}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-download

gathering-databases-download-modified: ${DATABASES_DOWNLOAD}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-download-modified

gathering-databases-scrape: ${DATABASES_SCRAPE}
	mkdir -p ${INTERIM_DB_PATH}
	make -C ${SRC_GATHERING_DB_PATH} gathering-databases-scrape

gathering-translation-quick: gathering-translation-common-quick gathering-translation-tcm-quick

gathering-translation-full: gathering-pmcid gathering-gbif gathering-chinese-board gathering-translation-tcmid gathering-translation-common gathering-translation-tcm

gathering-translation-common:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-common

gathering-translation-common-quick:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-common-quick

gathering-translation-tcm:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-tcm

gathering-translation-tcm-quick:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-tcm-quick

gathering-translation-tcmid:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-translation-tcmid

gathering-gbif:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-gbif

gathering-chinese-board:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-chinese-board

gathering-pmcid:
	make -C ${SRC_GATHERING_TRANSLATION_PATH} gathering-pmcid

gathering-taxonomy-full: gathering-taxonomy-npclassifier gathering-taxonomy-otl gathering-taxonomy-classyfire

gathering-taxonomy-quick: gathering-taxonomy-npclassifier gathering-taxonomy-classyfire

gathering-taxonomy-classyfire:
	make -C ${SRC_GATHERING_TAXONOMY_PATH} gathering-taxonomy-classyfire

gathering-taxonomy-npclassifier:
	make -C ${SRC_GATHERING_TAXONOMY_PATH} gathering-taxonomy-npclassifier

gathering-taxonomy-otl:
	make -C ${SRC_GATHERING_TAXONOMY_PATH} gathering-taxonomy-otl

curating-and-analyzing-and-visualizing: curating analyzing visualizing

curating-and-analyzing: curating analyzing

curating: get-bins curating-1-integrating curating-editing curating-3-integrating

curating-1-integrating: ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz: ${DATABASES} paths.mk ${SRC_PATH}/paths.R ${SRC_CURATING_PATH}/1_integrating.R
	cd src && Rscript ${SRC_CURATING_PATH}/1_integrating.R

curating-editing: curating-editing-structure curating-editing-organism curating-editing-reference

curating-editing-structure: curating-editing-structure-translating curating-editing-structure-integrating curating-editing-structure-sanitizing curating-editing-structure-stereocounting curating-editing-structure-naming curating-editing-structure-classifying

curating-editing-structure-translating: curating-editing-structure-translating-name curating-editing-structure-translating-inchi

curating-editing-structure-translating-name: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/names.R ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/nominal.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/names.R

curating-editing-structure-translating-inchi: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/inchi.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/inchi.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/inchi.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/inchi.tsv.gz
	cd src && python ${SRC_CURATING_EDITING_STRUCTURE_TRANSLATING_PATH}/inchi.py ${INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH}/inchi.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/inchi.tsv.gz structureOriginal_inchi

curating-editing-structure-integrating: ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz
${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/unique.tsv.gz: ${SRC_CURATING_EDITING_STRUCTURE_PATH}/2_integrating.R ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/inchi.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz
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
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz: $(wildcard ${INTERIM_TABLE_ORIGINAL_ORGANISM_PATH}/^[0-9]{6}.tsv) ${EXTERNAL_DICTIONARY_SOURCE_PATH}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_PATH}/1_processingOriginal.R
	cd src && Rscript ${SRC_CURATING_EDITING_ORGANISM_PATH}/1_processingOriginal.R

processing-organism-interim:
	-rm edit ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz

curating-editing-organism-translating: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz
${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build/libs/shadow.jar: ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build.gradle.kts $(wildcard ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/src/main/kotlin/*.kt)
	./gradlew castShadows
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/interim.tsv.gz: ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/deny.tsv ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/manualSubtraction.tsv ${SRC_CURATING_EDITING_ORGANISM_PATH}/2_translating_organism_kotlin/build/libs/shadow.jar #  ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/names.tsv.gz ${EXTERNAL_DICTIONARY_SOURCE_PATH}/tcm/names.tsv.gz 
	@java -jar ${SRC_CURATING_EDITING_ORGANISM_PATH_KT}/build/libs/shadow.jar ${DATA_PATH} ${MODE}

curating-editing-organism-processing-translated: curating-editing-organism-translating ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated.tsv.gz
${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated.tsv.gz: $(wildcard ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/translated/*.json) ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/original.tsv.gz ${EXTERNAL_DICTIONARY_SOURCE_PATH}/taxa/ranks.tsv ${SRC_CURATING_EDITING_ORGANISM_PATH}/3_processingTranslated.R
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
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/split.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/split.R $(wildcard ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/split/*.tsv.gz)
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/split.R

curating-editing-reference-translating-publishingDetails: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/publishingDetails.R $(wildcard ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/publishingDetails/*.tsv.gz)
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/publishingDetails.R

curating-editing-reference-translating-original: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/original.R $(wildcard ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/original/*.tsv.gz)
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_TRANSLATING_PATH}/original.R

curating-editing-reference-integrating: ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz
${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_PATH}/2_integrating.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/doi.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/pubmed.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/title.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/split.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/original.tsv.gz ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/publishingDetails.tsv.gz ${INTERIM_TABLE_ORIGINAL_REFERENCE_PATH}/full.tsv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_PATH}/2_integrating.R

curating-editing-reference-processing: ${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz
${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz: ${SRC_CURATING_EDITING_REFERENCE_PATH}/3_processing.R ${INTERIM_TABLE_TRANSLATED_REFERENCE_PATH}/integrated.tsv.gz # ${EXTERNAL_TRANSLATION_SOURCE_PATH}/pubmed/PMC-ids.csv.gz
	cd src && Rscript ${SRC_CURATING_EDITING_REFERENCE_PATH}/3_processing.R

curating-3-integrating: ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz: ${SRC_CURATING_PATH}/3_integrating.R ${INTERIM_TABLE_ORIGINAL_PATH}/table.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/inchi.tsv.gz ${INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH}/nominal.tsv.gz ${INTERIM_TABLE_PROCESSED_STRUCTURE_PATH}/processed.tsv.gz ${INTERIM_TABLE_PROCESSED_ORGANISM_PATH}/processed.tsv.gz ${INTERIM_TABLE_PROCESSED_REFERENCE_PATH}/processed.tsv.gz
	cd src && Rscript ${SRC_CURATING_PATH}/3_integrating.R

analyzing: gathering-validation	analyzing-sampling analyzing-validating analyzing-metrics analyzing-examples

analyzing-quick: gathering-validation analyzing-validating

analyzing-sampling: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/1_sampling.R

analyzing-validating: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/2_validating.R

analyzing-metrics: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz
	cd src && Rscript ${SRC_ANALYSING_PATH}/3_metrics.R

analyzing-examples: # ${INTERIM_TABLE_CURATED_PATH}/table.tsv.gz ${SRC_ANALYSING_PATH}/examples.R
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

lotus-quick: gathering-quick curating

lotus: gathering-full curating analyzing visualizing

lotus-full: gathering-full-hard curating analyzing visualizing

manual-entry: gathering-custom-dictionaries gathering-translation-quick gathering-taxonomy-quick curating analyzing-quick
