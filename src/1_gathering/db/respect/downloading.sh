#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/respect/respect ]; then

mkdir -p ../data/external/dbSource/respect/respect

wget "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip" -O ../data/external/dbSource/respect/respect.zip

unzip ../data/external/dbSource/respect/respect.zip -d ../data/external/dbSource/respect/respect/

fi

echo "Done"