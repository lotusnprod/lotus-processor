#!/usr/bin/env bash
# -*- coding: utf-8 -*-

eval $(../config.mk)

if [ ! -f ../data/external/taxonomySource/structure/npclasssifier/index_v"$INDEX_VERSION".json ]; then
  mkdir -p ../data/external/taxonomySource/structure/npclassifier/
  echo "Downloading"
  curl -o ../data/external/taxonomySource/structure/npclassifier/index_v"$INDEX_VERSION".json https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v"$INDEX_VERSION".json
fi

echo "You can now import in sql: "
echo " it is a json"
echo " next steps to do"
