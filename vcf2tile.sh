#!/bin/bash

source /opt/gcc-4.9.1/setup.sh
export PATH=/opt/python-2.7/bin:/cluster_share/GenomicsDB/GenomicsDB/bin:/cluster_share/GenomicsDB/samtools-1.3.1/bin:/cluster_share/GenomicsDB/bcftools-1.3.1/bin:$PATH
export LD_LIBRARY_PATH=/opt/gcc-4.9.1/lib64:/opt/gcc-4.9.1/lib:/usr/lib64/openmpi/lib
export HOME=/cluster_share/GenomicsDB/

for i in "$@"
do
case $i in
  -i=*|--inputs=*)
  INPUTS="${i#*=}"

  ;;
  
  -s|--sampleTag)
  SAMPLETAG=-s

  ;;

  -L|--load)
  LOAD=-L

  ;;

  -t|--tar)
  TAR=-t
  ;;
  *)

  ;;

esac
done

python /cluster_share/GenomicsDB/vcf2tile/vcf2tile.py -l /cluster_share/GenomicsDB/vcf2tile/loader.config -i ${INPUTS} $SAMPLETAG $LOAD $TAR

