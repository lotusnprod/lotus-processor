include	../../../paths.mk

.PHONY:	gathering-translation-common	gathering-translation-tcm
.PRECIOUS: %.tsv %.zip	%.json	%.gz

gathering-translation-common:	${INTERIM_DICTIONARY_PATH}/common/names.tsv.gz
${INTERIM_DICTIONARY_PATH}/common/names.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv	${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv	${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv	${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv	${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current/Taxon.tsv.zip	${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current/VernacularName.tsv.zip	${SRC_GATHERING_TRANSLATION_PATH}/common.R ${INTERIM_DICTIONARY_PATH}/common/manualSubtraction.tsv
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_TRANSLATION_PATH}/common.R

gathering-translation-tcm:	${INTERIM_DICTIONARY_PATH}/tcm/names.tsv.gz
${INTERIM_DICTIONARY_PATH}/tcm/names.tsv.gz:	${INTERIM_DICTIONARY_PATH}/tcm/manualSubtraction.tsv	${INTERIM_DICTIONARY_PATH}/latin/plantParts.tsv	${INTERIM_DICTIONARY_PATH}/latin/genitive/is.tsv	${INTERIM_DICTIONARY_PATH}/latin/genitive/i.tsv	${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/Chinese-Medicine-Board---List---Nomenclature-list-of-commonly-used-Chinese-herbal-medicines.XLSX	${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/tcmid/data/herb-TCMID.v2.01.txt	${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx	${SRC_GATHERING_TRANSLATION_PATH}/tcm.R
	cd	${PWD}/src	&&	Rscript	${SRC_GATHERING_TRANSLATION_PATH}/tcm.R