#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/afrotryp/44_2016_1764_MOESM1_ESM.doc ]; then

mkdir -p ../data/external/dbSource/afrotryp/

wget -nv "https://static-content.springer.com/esm/art%3A10.1007%2Fs00044-016-1764-y/MediaObjects/44_2016_1764_MOESM1_ESM.doc" -O ../data/external/dbSource/afrotryp/44_2016_1764_MOESM1_ESM.doc

fi

echo "Done"