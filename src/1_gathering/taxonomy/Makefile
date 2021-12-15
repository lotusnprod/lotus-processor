.PHONY: gathering-taxonomy-classyfire gathering-taxonomy-npclassifier gathering-taxonomy-otl
.PRECIOUS: %.tsv %.json %.tgz

gathering-taxonomy-classyfire: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/classyfire/tax_nodes.json
${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/classyfire/tax_nodes.json:
	@[ ! -f ${PWD}/data/external/taxonomySource/structure/classyfire/tax_nodes.json ] && mkdir -p ${PWD}/data/external/taxonomySource/structure/classyfire/ && echo "Downloading" && curl -o ${PWD}/data/external/taxonomySource/structure/classyfire/tax_nodes.json http://classyfire.wishartlab.com/tax_nodes.json || true

gathering-taxonomy-npclassifier: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/npclassfifier/index_v${INDEX_VERSION}.json
${EXTERNAL_TAXONOMY_SOURCE_PATH}/structure/npclassfifier/index_v${INDEX_VERSION}.json:
	@[ ! -f ${PWD}/data/external/taxonomySource/structure/npclassifier/index_v${INDEX_VERSION}.json ] && mkdir -p ${PWD}/data/external/taxonomySource/structure/npclassifier/ && echo "Downloading" && curl -o ${PWD}/data/external/taxonomySource/structure/npclassifier/index_v${INDEX_VERSION}.json https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v${INDEX_VERSION}.json || true

gathering-taxonomy-otl: ${EXTERNAL_TAXONOMY_SOURCE_PATH}/organism/taxonomy.tsv
${EXTERNAL_TAXONOMY_SOURCE_PATH}/organism/taxonomy.tsv:
	@[ ! -f ${PWD}/data/external/taxonomySource/organism/ott${OTT_VERSION}.tgz ] && mkdir -p ${PWD}/data/external/taxonomySource/organism/ && echo "Downloading" && curl -o ${PWD}/data/external/taxonomySource/organism/ott${OTT_VERSION}.tgz https://files.opentreeoflife.org/ott/ott${OTT_VERSION}/ott${OTT_VERSION}.tgz && echo "Extracting" && tar -xzf ${PWD}/data/external/taxonomySource/organism/ott${OTT_VERSION}.tgz --strip-components=1 ott${OTT_VERSION}/taxonomy.tsv && sed "s/\t|\t/\t/g" <${PWD}/taxonomy.tsv | sed "s/\t$$//g" >${PWD}/data/external/taxonomySource/organism/taxonomy.tsv && rm ${PWD}/taxonomy.tsv  && echo "You can now import in sqlite: " && echo " create if not exist table ott_taxonomy(uid int, parent_uid int, name text, rank text, sourceinfo text, uniqname text, flags text);" && echo " .separator \"\\t\"" && echo " .import ott_taxonomy.tsv ott_taxonomy" || true