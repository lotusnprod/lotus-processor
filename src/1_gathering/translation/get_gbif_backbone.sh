#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/translationSource/common/backbone-current.zip ]; then

mkdir -p ../data/external/translationSource/common/

# wget -nv "https://hosted-datasets.gbif.org/datasets/backbone/backbone-current.zip" -O ../data/external/translationSource/common/backbone-current.zip

wget -nv "https://hosted-datasets.gbif.org/datasets/backbone/2021-03-03/backbone.zip" -O ../data/external/translationSource/common/backbone-current.zip

fi

echo "Done"