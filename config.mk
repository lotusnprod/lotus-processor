# To run in full mode for the build turn that to full
# It can be done at run time by doing
# make MODE=full <target>
export MODE ?= min

export DATA_PATH ?= ${PWD}/data
export SRC_PATH ?= ${PWD}/src
export BIN_PATH ?= ${PWD}/bin
export TESTS_PATH ?= ${PWD}/tests

export GNFINDER_VERSION = v0.16.2
export GNVERIFIER_VERSION = v0.5.2
export OPSIN_VERSION = 2.5.0

export OTT_VERSION = 3.3

export NPCLASSIFIER_VERSION = 1.5
export INDEX_VERSION = 1

export UNAME := $(shell uname)

PLATFORM := unsupported
ifeq ($(OS),Windows_NT)
	PLATFORM := unsupported
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        PLATFORM := linux
    endif
    ifeq ($(UNAME_S),Darwin)
        PLATFORM := mac
    endif
endif
