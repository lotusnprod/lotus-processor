.PHONY: gathering-custom-dictionaries
.PRECIOUS: %.zip %.tsv %.gz

gathering-custom-dictionaries: ${EXTERNAL_DICTIONARY_SOURCE_PATH}/dbNames.tsv
${EXTERNAL_DICTIONARY_SOURCE_PATH}/dbNames.tsv:
	@[ ! -f ${EXTERNAL_DICTIONARY_SOURCE_PATH} ] && mkdir -p ${EXTERNAL_DICTIONARY_SOURCE_PATH} && echo "Downloading" && wget -nv "https://zenodo.org/record/${ZENODO_CUSTOM_DIC}/files/dictionarySource.zip?download=1" -O ${PWD}/data/external/dictionarySource.zip && unzip ${PWD}/data/external/dictionarySource.zip -d ${PWD}/data/external/ && rm ${PWD}/data/external/dictionarySource.zip && echo "Done" || true
