#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/coconut/COCONUT_DB.sdf ]; then

mkdir -p ../data/external/dbSource/coconut/

wget -nv "https://coconut.naturalproducts.net/download/sdf" -O ../data/external/dbSource/coconut/COCONUT_DB.sdf

fi

echo "Done"