#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dictionarySource/ ]; then

mkdir -p ../data/external/dictionarySource/

wget -nv "https://zenodo.org/record/5801816/files/dictionarySource.zip?download=1" -O ../data/external/dictionarySource.zip

unzip ../data/external/dictionarySource.zip -d ../data/external/

rm ../data/external/dictionarySource.zip

fi

echo "Done"