#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/dictionarySource/

wget "https://osf.io/gvre8/download" -O ../data/external/dictionarySource.zip

unzip ../data/external/dictionarySource.zip -d ../data/external/

rm ../data/external/dictionarySource.zip