#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/translationSource/pubmed/

wget "https://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz" -O ../data/external/translationSource/pubmed/PMC-ids.csv.gz
