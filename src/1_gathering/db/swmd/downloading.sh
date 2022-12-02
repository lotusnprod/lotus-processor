#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/swmd/Mol ]; then

mkdir -p ../data/external/dbSource/swmd

wget -nv "http://www.swmd.co.in/Mol.zip" -O ../data/external/dbSource/swmd/Mol.zip

unzip ../data/external/dbSource/swmd/Mol.zip -d ../data/external/dbSource/swmd/

fi

echo "Done"