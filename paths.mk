export DATA_PATH ?= ${PWD}/data
#export DATA_PATH ?= ${PWD}/data_min
export SRC_PATH ?= ${PWD}/src


export INTERIM_DB_PATH = ${DATA_PATH}/interim/db
export INTERIM_DICTIONARY_PATH = ${DATA_PATH}/interim/dictionaries
export DB_SOURCE_PATH = ${DATA_PATH}/external/dbSource
## be carfeul with comments they can ruin the makefile as stated here: https://www.gnu.org/software/make/manual/html_node/Recipe-Syntax.html. At least it made mine do not work anymore

# Below is and endless PATH nightmare. Why dont we have a unique ITERIM folder to dump all files, uniquely named ?
## This is because the files should be slowly removed and renamed in dictionaries instead. so working with interim/original table at the begining, bunch of dictionaries (original to translated, translated to curated etc) and processed table in the end, in my view at least
# export INTERIM_TABLE_PATH = ${DATA_PATH}/interim/tables
export INTERIM_TABLE_PATH = ${DATA_PATH}/interim/tables_min
export INTERIM_TABLE_ORIGINAL_PATH = ${INTERIM_TABLE_PATH}/0_original
export INTERIM_TABLE_ORIGINAL_ORGANISM_PATH = ${INTERIM_TABLE_ORIGINAL_PATH}/organism
export INTERIM_TABLE_ORIGINAL_REFERENCE_PATH = ${INTERIM_TABLE_ORIGINAL_PATH}/reference
export INTERIM_TABLE_ORIGINAL_STRUCTURE_PATH = ${INTERIM_TABLE_ORIGINAL_PATH}/structure
export INTERIM_TABLE_TRANSLATED_PATH = ${INTERIM_TABLE_PATH}/1_translated
export INTERIM_TABLE_TRANSLATED_STRUCTURE_PATH = ${INTERIM_TABLE_TRANSLATED_PATH}/structure
export INTERIM_TABLE_CLEANED_PATH = ${INTERIM_TABLE_PATH}/2_cleaned
export INTERIM_TABLE_CLEANED_ORGANISM_PATH = ${INTERIM_TABLE_CLEANED_PATH}/organism
export INTERIM_TABLE_CLEANED_REFERENCE_PATH = ${INTERIM_TABLE_CLEANED_PATH}/reference
export INTERIM_TABLE_CLEANED_STRUCTURE_PATH = ${INTERIM_TABLE_CLEANED_PATH}/structure

export SRC_GATHERING_PATH = ${SRC_PATH}/1_gathering/db
export SRC_CURATING_PATH = ${SRC_PATH}/2_curating

export SRC_CURATING_1_INTEGRATING_PATH = ${SRC_CURATING_PATH}/1_integrating
export SRC_CURATING_EDITING_PATH = ${SRC_CURATING_PATH}/2_editing
export SRC_CURATING_EDITING_STRUCTURE_PATH =  ${SRC_CURATING_EDITING_PATH}/structure
export SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_PATH =  ${SRC_CURATING_EDITING_STRUCTURE_PATH}/subscripts
export SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_TRANSLATING_PATH =  ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_PATH}/1_translating
export SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_INTEGRATING_PATH =  ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_PATH}/2_integrating
export SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_CLEANINGANDENRICHING_PATH =  ${SRC_CURATING_EDITING_STRUCTURE_SUBSCRIPTS_PATH}/3_cleaningAndEnriching
export SRC_CURATING_3_INTEGRATING_PATH = ${SRC_CURATING_PATH}/3_integrating
