export DATA_PATH ?= ${PWD}/data
#export DATA_PATH ?= ${PWD}/data_min
export SRC_PATH ?= ${PWD}/src

export INTERIM_PATH = ${DATA_PATH}/interim/db
export DB_SOURCE_PATH = ${DATA_PATH}/external/dbSource
## be carfeul with comments they can ruin the makefile as stated here: https://www.gnu.org/software/make/manual/html_node/Recipe-Syntax.html. At least it made mine do not work anymore

# Below is and endless PATH nightmare. Why dont we have a unique ITERIM folder to dump all files, uniquely named ?
## This is because the files should be slowly removed and renamed in dictionaries instead. so working with interim/original table at the begining, bunch of dictionaries (original to translated, translated to curated etc) and processed table in the end, in my view at least
# export INTERIM_TABLE_PATH = ${DATA_PATH}/interim/tables
export INTERIM_TABLE_PATH = ${DATA_PATH}/interim/tables_min
export INTERIM_TABLE_ORIGINAL_PATH = ${INTERIM_TABLE_PATH}/0_original
export INTERIM_TABLE_TRANSLATED_PATH = ${INTERIM_TABLE_PATH}/1_translated
export INTERIM_TABLE_CLEANED_PATH = ${INTERIM_TABLE_PATH}/2_cleaned
export INTERIM_TABLE_CURATED_PATH = ${INTERIM_TABLE_PATH}/3_curated

export SRC_GATHERING_PATH = ${SRC_PATH}/1_gathering/db
export SRC_CURATING_PATH = ${SRC_PATH}/2_curating

export SRC_CURATING_EDITING_PATH = ${SRC_CURATING_PATH}/2_editing
export SRC_CURATING_EDITING_CHEMO_PATH =  ${SRC_CURATING_EDITING_PATH}/chemo
export SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_PATH =  ${SRC_CURATING_EDITING_CHEMO_PATH}/subscripts
export SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_TRANSLATING_PATH =  ${SRC_CURATING_EDITING_CHEMO_SUBSCRIPTS_PATH}/1_translating
