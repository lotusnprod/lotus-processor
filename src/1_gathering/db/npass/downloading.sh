#!/usr/bin/env bash
# -*- coding: utf-8 -*-

if [ ! -d ../data/external/dbSource/npass/ ]; then

mkdir -p ../data/external/dbSource/npass/

wget -nv "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_generalInfo.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt

wget -nv "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_properties.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_properties.txt

wget -nv "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_speciesInfo.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt

wget -nv "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_species_pair.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_species_pair.txt

fi

echo "Done"