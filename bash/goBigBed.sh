#!/usr/bin/env bash
#################################################################
#
#date: June 22, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# goBigBed.sh
# Usage : 
# goBigBed.sh -d Absolute Path to bedfiles
# Purpose :
# Transform bed to bigBed
#
# Note : 
# Path to GRCh38.primary_assembly.genome.sizes hard coded. Need to modify manually inside code.
#
#
#################################################################

while getopts d: option
do
        case "${option}"
        in
              d) DIR=${OPTARG};;

        esac
done
echo ${DIR}

echo "Sort"

find ${DIR} -name "*.bed" |  xargs -I{} bash -c 'sort -k1,1 -k2,2n $0 > ${0/.bed/.sorted.bed}' {} \;

echo "Bigbed"
find ${DIR} -name "*.sorted.bed" |  xargs -I{} bash -c 'bedToBigBed  -type=bed6 $0 /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genome/GRCh38_PRIM_GENCODE_R25/GRCh38.primary_assembly.genome.sizes ${0/.sorted.bed/.sorted.bb}' {} \;

ls ${DIR}
