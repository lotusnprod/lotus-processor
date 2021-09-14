.PHONY:	gathering-translation-common	gathering-translation-tcm
.PRECIOUS: %.tsv %.zip	%.json	%.gz

gathering-translation-common:	${INTERIM_DICTIONARY_PATH_FIX}/common/names.tsv.gz
${INTERIM_DICTIONARY_PATH_FIX}/common/names.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv	${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv	${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv	${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv	${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current.zip		${SRC_GATHERING_TRANSLATION_PATH}/common.R ${INTERIM_DICTIONARY_PATH_FIX}/common/manualSubtraction.tsv
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_TRANSLATION_PATH}/common.R

gathering-translation-tcm:	${INTERIM_DICTIONARY_PATH_FIX}/tcm/names.tsv.gz
${INTERIM_DICTIONARY_PATH_FIX}/tcm/names.tsv.gz:	${INTERIM_DICTIONARY_PATH_FIX}/tcm/manualSubtraction.tsv	${INTERIM_DICTIONARY_PATH_FIX}/latin/plantParts.tsv	${INTERIM_DICTIONARY_PATH_FIX}/latin/genitive/is.tsv	${INTERIM_DICTIONARY_PATH_FIX}/latin/genitive/i.tsv	${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/Chinese-Medicine-Board---List---Nomenclature-list-of-commonly-used-Chinese-herbal-medicines.XLSX	${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/tcmid/data/herb-TCMID.v2.01.txt	${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx	${SRC_GATHERING_TRANSLATION_PATH}/tcm.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_TRANSLATION_PATH}/tcm.R