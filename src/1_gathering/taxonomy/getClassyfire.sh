#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/taxonomySource/structure/classyfire/tax_nodes.json ]; then
  mkdir -p ../data/external/taxonomySource/structure/classyfire/
  echo "Downloading"
  curl -o ../data/external/taxonomySource/structure/classyfire/tax_nodes.json http://classyfire.wishartlab.com/tax_nodes.json
fi

echo "You can now import in sql: "
echo " it is a json"
echo " next steps to do"