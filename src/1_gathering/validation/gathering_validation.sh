#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/validation/ ]; then

mkdir -p ../data/validation/dictionarySource/

wget -nv "https://osf.io/vg2we/download" -O ../data/validation.zip

unzip ../data/validation.zip -d ../data/

rm ../data/validation.zip

fi

echo "Done"