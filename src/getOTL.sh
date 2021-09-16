#!/usr/bin/env bash
# -*- coding: utf-8 -*-

OTT_VERSION="3.3"

if [ ! -f ../data/external/taxonomySource/organism/ott"$OTT_VERSION".tgz ]; then
  echo "Downloading"
  curl -o ../data/external/taxonomySource/organism/ott"$OTT_VERSION".tgz https://files.opentreeoflife.org/ott/ott"$OTT_VERSION"/ott"$OTT_VERSION".tgz
fi

if [ ! -f ott_taxonomy.tsv ]; then
  echo "Extracting"
  tar -xzf ../data/external/taxonomySource/organism/ott"$OTT_VERSION".tgz --strip-components=1 ott"$OTT_VERSION"/taxonomy.tsv
  sed "s/\t|\t/\t/g" <taxonomy.tsv | sed "s/\t$//g" >../data/external/taxonomySource/organism/taxonomy.tsv
  rm taxonomy.tsv
fi

echo "You can now import in sqlite: "
echo " create if not exist table ott_taxonomy(uid int, parent_uid int, name text, rank text, sourceinfo text, uniqname text, flags text);"
echo " .separator \"\\t\""
echo " .import ott_taxonomy.tsv ott_taxonomy"
