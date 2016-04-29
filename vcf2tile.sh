#!/bin/bash

source /opt/gcc-4.9.1/setup.sh
export PATH=/opt/python-2.7/bin:/mnt/app_hdd/TileDB/GenomicsDB/bin:$PATH

python vcf2tile.py -l $1 -i $2 -s -L
