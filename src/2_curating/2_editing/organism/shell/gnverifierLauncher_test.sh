#!/bin/sh

# A simple script to use gnverifier

source 2_curating/2_editing/organism/shell/gnverifier.sh

gnverifier ../data/interim/tables_test/2_cleaned/organism/verify.tsv.gz ../data/interim/tables_test/2_cleaned/organism/verified.json
