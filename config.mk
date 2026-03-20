# To run in full mode for the build turn that to full
# It can be done at run time by doing
# make MODE=full <target>
export MODE ?= min

export DATA_PATH ?= ${PWD}/data
export SRC_PATH ?= ${PWD}/src
export BIN_PATH ?= ${PWD}/bin
export TESTS_PATH ?= ${PWD}/tests

export GNFINDER_VERSION = v1.1.6
export GNVERIFIER_VERSION = v1.3.1
export OPSIN_VERSION = 2.9.0

export OTT_VERSION = 3.7.3
export GBIF_BACKBONE = 2023-08-28

export NPCLASSIFIER_VERSION = 1.5
export INDEX_VERSION = 1

export OSF_VALIDATION = vg2we
export ZENODO_CUSTOM_DIC = 6487114

export UNAME := $(shell uname)

# get OS info
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# normalize OS
ifeq ($(UNAME_S),Windows_NT)
  OS := win
else ifeq ($(UNAME_S),Linux)
  OS := linux
else ifeq ($(UNAME_S),Darwin)
  OS := mac
else
  OS := unsupported
endif

# normalize architecture
ifeq ($(UNAME_M),x86_64)
  ARCH := amd64
else ifeq ($(UNAME_M),aarch64)
  ARCH := arm64
else ifeq ($(UNAME_M),arm64)
  ARCH := arm64
else
  ARCH := unknown
endif

PLATFORM := $(OS)-$(ARCH)

NPROCS := 1
ifeq ($(UNAME_S),Linux)
  NPROCS := $(shell grep -c ^processor /proc/cpuinfo)
endif
ifeq ($(UNAME_S),Darwin)
  NPROCS := $(shell sysctl -n hw.ncpu)
endif

show-platform:
	@echo $(PLATFORM)
