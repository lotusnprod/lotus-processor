#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/cyanometdb/CyanoMetDB_WR_Feb2021.csv ]; then

mkdir -p ../data/external/dbSource/cyanometdb/

## Thanks to the Norman suspect list
wget "https://zenodo.org/record/4562688/files/CyanoMetDB_WR_Feb2021.csv?download=1" -O ../data/external/dbSource/cyanometdb/CyanoMetDB_WR_Feb2021.csv

fi

echo "Done"