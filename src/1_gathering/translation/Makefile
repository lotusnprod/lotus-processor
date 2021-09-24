.PHONY: gathering-translation-common gathering-translation-tcm gathering-gbif gathering-chinese-board gathering-pmcid gathering-translation-tcmid
.PRECIOUS: %.csv %.tsv %.zip %.json %.gz %.XLSX %.txt %.xlsx

gathering-translation-common: ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/names.tsv.gz
${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/names.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv ${EXTERNAL_DB_SOURCE_PATH}/foodb/foodb_2020_04_07_csv/Food_copy.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current.zip ${SRC_GATHERING_TRANSLATION_PATH}/common.R ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/manualSubtraction.tsv
	cd ${PWD}/src && Rscript ${SRC_GATHERING_TRANSLATION_PATH}/common.R

gathering-translation-common-quick: ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/names.tsv.gz
${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/names.tsv.gz: ${EXTERNAL_DB_SOURCE_PATH}/phenolexplorer/foods.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/COMMON_NAMES.csv ${EXTERNAL_DB_SOURCE_PATH}/drduke/Duke-Source-CSV/FNFTAX.csv ${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current.zip ${SRC_GATHERING_TRANSLATION_PATH}/common.R ${EXTERNAL_DICTIONARY_SOURCE_PATH}/common/manualSubtraction.tsv
	cd ${PWD}/src && Rscript ${SRC_GATHERING_TRANSLATION_PATH}/common.R

gathering-translation-tcm: ${EXTERNAL_DICTIONARY_SOURCE_PATH}/tcm/names.tsv.gz
${EXTERNAL_DICTIONARY_SOURCE_PATH}/tcm/names.tsv.gz: ${EXTERNAL_DICTIONARY_SOURCE_PATH}/tcm/manualSubtraction.tsv ${EXTERNAL_DICTIONARY_SOURCE_PATH}/latin/plantParts.tsv ${EXTERNAL_DICTIONARY_SOURCE_PATH}/latin/genitive/is.tsv ${EXTERNAL_DICTIONARY_SOURCE_PATH}/latin/genitive/i.tsv ${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX ${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/tcmid/data/herb-TCMID.v2.01.txt ${EXTERNAL_DB_SOURCE_PATH}/tmmc/compound.xlsx ${SRC_GATHERING_TRANSLATION_PATH}/tcm.R
	cd ${PWD}/src && Rscript ${SRC_GATHERING_TRANSLATION_PATH}/tcm.R

gathering-gbif: ${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current.zip
${EXTERNAL_TRANSLATION_SOURCE_PATH}/common/backbone-current.zip: ${SRC_GATHERING_TRANSLATION_PATH}/get_gbif_backbone.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TRANSLATION_PATH}/get_gbif_backbone.sh

gathering-chinese-board: ${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX
${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/Chinese-Medicine-Board---List---Nomenclature-compendium-of-commonly-used-Chinese-herbal-medicines.XLSX: ${SRC_GATHERING_TRANSLATION_PATH}/get_chinese_medicine_board.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TRANSLATION_PATH}/get_chinese_medicine_board.sh

gathering-pmcid: ${EXTERNAL_TRANSLATION_SOURCE_PATH}/pubmed/PMC-ids.csv.zip
${EXTERNAL_TRANSLATION_SOURCE_PATH}/pubmed/PMC-ids.csv.zip: ${SRC_GATHERING_TRANSLATION_PATH}/get_pmcid.sh
	cd ${PWD}/src && bash ${SRC_GATHERING_TRANSLATION_PATH}/get_pmcid.sh

gathering-translation-tcmid: # $(wildcard ${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/tcmid/**)
# $(wildcard ${EXTERNAL_TRANSLATION_SOURCE_PATH}/tcm/tcmid/**):	
	cd ${PWD}/src && bash ${SRC_GATHERING_TRANSLATION_PATH}/get_tcmid.sh
