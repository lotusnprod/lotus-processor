#!/bin/sh

# A simple script to use gnverifier

source 2_curating/2_editing/organism/shell/gnverifier_capital.sh

gnverifier_capital ../data/interim/tables/0_original/organism/original.tsv.gz ../data/interim/tables/2_cleaned/organism/original_verified.json
