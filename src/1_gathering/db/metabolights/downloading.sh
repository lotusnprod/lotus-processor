#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -f ../data/external/dbSource/metabolights/eb-eye_metabolights_complete.xml ]; then

mkdir -p ../data/external/dbSource/metabolights/

wget -r "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb-eye_metabolights_complete.xml" -O ../data/external/dbSource/metabolights/eb-eye_metabolights_complete.xml

rm -r "ftp.ebi.ac.uk"

fi

echo "Done"