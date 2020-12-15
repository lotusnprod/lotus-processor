# To run in full mode for the build turn that to full
# It can be done at run time by doing
# make MODE=full <target>
export MODE ?= min

export DATA_PATH ?= ${PWD}/data
export SRC_PATH ?= ${PWD}/src
export BIN_PATH ?= ${PWD}/bin