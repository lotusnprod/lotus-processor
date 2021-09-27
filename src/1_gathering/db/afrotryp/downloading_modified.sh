#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/afrotryp/afrotryp.tsv.zip ]; then

	mkdir -p ../data/external/dbSource/afrotryp/

	wget -nv "https://osf.io/q47fd/download" -O ../data/external/dbSource/afrotryp/afrotryp.tsv.zip
fi

echo "Done"
