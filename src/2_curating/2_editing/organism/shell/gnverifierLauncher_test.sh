#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# A simple script to use gnverifier

source 2_curating/2_editing/organism/shell/gnverifier.sh

gnverifier ../data/interim/tables_test/2_processed/organism/verify.tsv.gz ../data/interim/tables_test/2_processed/organism/verified.json
