#!/usr/bin/env bash
# -*- coding: utf-8 -*-

mkdir -p ../data/external/dbSource/npass/

wget "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_generalInfo.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_generalInfo.txt

wget "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_properties.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_properties.txt

wget "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_speciesInfo.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_speciesInfo.txt

wget "http://bidd.group/NPASS/downloadFiles/NPASSv1.0_download_naturalProducts_species_pair.txt" -O ../data/external/dbSource/npass/NPASSv1.0_download_naturalProducts_species_pair.txt

