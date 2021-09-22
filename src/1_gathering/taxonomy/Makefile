.PHONY: gathering-taxonomy-classyfire gathering-taxonomy-npclassifier gathering-taxonomy-otl
.PRECIOUS: %.tsv %.json %.tgz

gathering-taxonomy-classyfire: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/classyfire/tax_nodes.json
${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/classyfire/tax_nodes.json: ${SRC_GATHERING_TAXONOMY_PATH}/getClassyfire.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TAXONOMY_PATH}/getClassyfire.sh

gathering-taxonomy-npclassifier: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/npclassfifier/index_v${INDEX_VERSION}.json
${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/npclassfifier/index_v${INDEX_VERSION}.json: ${SRC_GATHERING_TAXONOMY_PATH}/getNPClassifier.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TAXONOMY_PATH}/getNPClassifier.sh

gathering-taxonomy-otl: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/organism/taxonomy.tsv
${EXTERNAL_TAXONOMY_SOURCE_PATH}/organism/taxonomy.tsv: ${SRC_GATHERING_TAXONOMY_PATH}/getOTL.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TAXONOMY_PATH}/getOTL.sh
