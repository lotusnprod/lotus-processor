#!/usr/bin/env bash

echo "Installing packages missing from conda"

conda skeleton cran tidyjson --output-dir $CONDA/envs/lotus_env/lib/R/library
conda skeleton cran classyfirer --allow-archived --output-dir $CONDA/envs/lotus_env/lib/R/library
conda skeleton cran https://github.com/ropensci/rcrossref --output-dir $CONDA/envs/lotus_env/lib/R/library