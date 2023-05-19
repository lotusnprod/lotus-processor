#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/cyanometdb/CyanoMetDB_v02_2023.csv ]; then

mkdir -p ../data/external/dbSource/cyanometdb/

## Thanks to the Norman suspect list
wget -nv "https://zenodo.org/record/7922070/files/CyanoMetDB_v02_2023.csv?download=1" -O ../data/external/dbSource/cyanometdb/CyanoMetDB_v02_2023.csv

fi

echo "Done"