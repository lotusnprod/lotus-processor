#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/npatlas/NPAtlas_download.tsv ]; then

mkdir -p ../data/external/dbSource/npatlas/

wget "https://www.npatlas.org/static/downloads/NPAtlas_download.tsv" -O ../data/external/dbSource/npatlas/NPAtlas_download.tsv

fi

echo "Done"