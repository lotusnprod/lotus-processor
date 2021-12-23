#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# A simple script to use gnverifier

source 2_curating/2_editing/organism/shell/gnverifier_capital.sh

gnverifier_capital ../data/interim/tables_custom/0_original/organism/original.tsv.gz ../data/interim/tables_custom/2_processed/organism/original_verified.json
