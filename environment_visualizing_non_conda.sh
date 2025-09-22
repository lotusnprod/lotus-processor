#!/usr/bin/env bash

echo "Installing packages missing from conda"

conda skeleton cran ggstar philentropy rotl rpostgresql
