#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/carotenoiddb/Carotenoids_InChI_InChIKey.tsv ]; then

mkdir -p ../data/external/dbSource/carotenoiddb/

wget -nv "http://carotenoiddb.jp/FTP/Carotenoids_InChI_InChIKey.tsv" -O ../data/external/dbSource/carotenoiddb/Carotenoids_InChI_InChIKey.tsv

fi

echo "Done"