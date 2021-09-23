#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/datawarrior/NaturalProductsNaturalProducts.txt ]; then

mkdir -p ../data/external/dbSource/datawarrior/

wget "https://osf.io/3g75a/download" -O ../data/external/dbSource/datawarrior/NaturalProductsNaturalProducts.txt

fi

echo "Done"