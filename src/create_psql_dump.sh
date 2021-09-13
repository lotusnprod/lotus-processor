#!/usr/bin/env bash
# -*- coding: utf-8 -*-

pg_dump -Fc -a -O -Z 9 lotus > ../data/processed/lotus.psql
