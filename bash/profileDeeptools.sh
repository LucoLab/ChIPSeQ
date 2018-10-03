#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"

#################################################################
#
#date: Jan 29, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# profileDeeptools.sh
# Usage : 
# profileDeeptools.sh & 12 parameters (see deeptools_countsignal)
# Purpose :
# This script is not used alone.
# It's called inside python script deeptools_countsignal.py
# It does deeptools stuff to create wig normalise and plot profile
#
#################################################################


echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
echo $7
echo $8
echo ${9}
echo ${10}
echo ${11}

rep1BamPath=$1
rep2BamPath=$2
controlBamPath=$3
rep1BamTotalReadMapped=$4
rep2BamTotalReadMapped=$5
controlBamTotalReadMapped=$6
outPath=$7
joinedPathToBed=$8
winInExon=$9
winInIntron=${10}
name=${11}
fragSize=${12}
duplicates=${13}

lengthFile=/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/length.chromosomes.sorted.bed
chromLengthIndex=/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/indexes/GRCh38_PRIM_GENCODE_R25/chrNameLength.txt
echo "WELCOME TO THE BASH"
echo ${outPath}
echo ${name}
echo ${fragSize}


OUT=/home/jean-philippe.villemin/data/data/bigWigDeeptools/

suffix=".bam";
rep1Path=${rep1BamPath%$suffix}; #Remove suffix
rep2Path=${rep2BamPath%$suffix}; #Remove suffix
controlPath=${controlBamPath%$suffix}; #Remove suffix

joinedPathToBed=${joinedPathToBed//[,]/ }

scaleRatio1=1000000/rep1BamTotalReadMapped
scaleRatio2=1000000/rep2BamTotalReadMapped 
scaleRatioControl=1000000/controlBamTotalReadMapped

if [ "${duplicates}" == "YES" ]; then
	duplicates=""
else
	duplicates="--ignoreDuplicates"
fi
# Default
duplicates=""

echo "DUP is ${duplicates}"
#fragSize=200


if [ -e ${OUT}${name}.rep1.new.bw ]
then
    echo "${OUT}${name}.rep1.new.bw Exists."
else	
    echo "Create bigwig rep1 " #--normalizeUsingRPKM --binSize 10 ignoreDuplicates True ? --smoothLength 30
bamCoverage --bam ${rep1BamPath}  -p 16 --outFileName ${OUT}${name}.rep1.new.bw ${duplicates} --normalizeUsingRPKM  --smoothLength 30 --binSize 10 --outFileFormat bigwig  --extendRead ${fragSize} --minMappingQuality 20
	#--smoothLength : raoul used it...
	#--scaleFactor I could set total number mapped /1000000 and remove normaliseUsingRPKM
fi

if [ -e ${OUT}${name}.rep2.new.bw ]
then
    echo "${OUT}${name}.rep2.new.bw Exists."
else
	echo "Create bigwig rep2 " #--normalizeUsingRPKM --binSize 10
bamCoverage --bam ${rep2BamPath}  -p 16  --outFileName ${OUT}${name}.rep2.new.bw ${duplicates} --normalizeUsingRPKM  --smoothLength 30 --binSize 10 --outFileFormat bigwig  --extendRead ${fragSize} --minMappingQuality 20
fi

if [ -e ${OUT}${name}.control.new.bw ]
then
    echo "${OUT}${name}.control.new.bw Exists."
else
	echo "Create bigwig control "  #--normalizeUsingRPKM
bamCoverage --bam ${controlBamPath} -p 16  --outFileName ${OUT}${name}.control.new.bw ${duplicates} --normalizeUsingRPKM  --smoothLength 30 --binSize 10 --outFileFormat bigwig  --extendRead ${fragSize} --minMappingQuality 20
fi

if [ -e ${OUT}${name}.mean.normalised.wig ]
then
    echo "${OUT}${name}.mean.normalised.wig Exists."
else
	echo "wiggletools mean Replicates and substract new"
	echo "wiggletools write ${OUT}${name}.mean.normalised.wig trim ${lengthFile} gt 0 diff mean ${OUT}${name}.rep1.new.bw ${OUT}${name}.rep2.new.bw : ${OUT}${name}.control.new.bw"
	wiggletools write ${OUT}${name}.mean.normalised.intermediate.wig  trim ${lengthFile} diff mean ${OUT}${name}.rep1.new.bw ${OUT}${name}.rep2.new.bw : ${OUT}${name}.control.new.bw
	/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if (  $4 < 0 )$4=0; print  $1,$2,$3,$4 ; }'  ${OUT}${name}.mean.normalised.intermediate.wig > ${OUT}${name}.mean.normalised.wig

	
fi

if [ -e ${OUT}${name}.mean.normalised.bw ]
then
    echo "${OUT}${name}.mean.normalised.bw Exists."
else
	echo "WigGoBigWig"
	echo "wigToBigWig -clip ${OUT}${name}.mean.normalised.wig ${chromLengthIndex} ${OUT}${name}.mean.normalised.bw"
	wigToBigWig -clip ${OUT}${name}.mean.normalised.wig ${chromLengthIndex} ${OUT}${name}.mean.normalised.bw
fi

if [ -e ${OUT}${name}.mean.normalised.gz ]
then
    echo "${OUT}${name}.mean.normalised.gz Exists."
else
	echo "computeMatrix"
	echo "	computeMatrix reference-point  -R ${joinedPathToBed} -S ${OUT}${name}.mean.normalised.bw --outFileName ${OUT}${name}.mean.normalised.gz -a ${winInExon} -b ${winInIntron}"
	computeMatrix reference-point  -R ${joinedPathToBed} --binSize 10 -S ${OUT}${name}.mean.normalised.bw --outFileName ${OUT}${name}.mean.normalised.gz -a ${winInExon} -b ${winInIntron}
	#-binSize ( Already binsize was set to 10 in bamCoverage , so no need to set it again,default is 10) 
	#-averageTypeBins (by default it use mean)
	#Length, in bases, of the non-overlapping bins for averaging the score over the regions length.
	#The bins used in computeMatrix will often only partially overlap those in 
	#the underlying bigWig files. The computeMatrix bins generally affect 
	#the resolution of the resulting image whereas those in the bigWig 
	#files more directly affect the underlying data. It doesn't make sense 
	#to use a bin size of 1 in computeMatrix when you've used a bin size of 
	#50 in bamCoverage/bamCompare, of course. The defaults (50 for bigWig 
	#files and 10 in computeMatrix) represent a convenient trade off, where 
	#space is saved in the bigWig files but you can still get some decent 
	#resolution in the resulting images. A bin averaging type of "mean" is 
	#the default for both. 
fi

if [ -e ${OUT}${name}.mean.normalised.png ]
then
	
    echo "${OUT}${name}.mean.normalised.png Exists."
else
	echo "plotProfile"
	echo "plotProfile -m ${OUT}${name}.mean.normalised.gz -out ${OUT}${name}.mean.normalised.png"
	
	plotProfile -m ${OUT}${name}.mean.normalised.gz -out ${OUT}${name}.mean.normalised.png
	
fi

rm ${OUT}${name}.rep1.new.bw
rm ${OUT}${name}.rep2.new.bw
rm ${OUT}${name}.control.new.bw
rm ${OUT}${name}.mean.normalised.wig 
rm ${OUT}${name}.mean.normalised.gz 
rm ${OUT}${name}.mean.normalised.intermediate.wig 

ARCHIVE_STORE=/home/jean-philippe.villemin/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_CHIPSEQ_FILES/normalized_bigwig/
EXTERN=/home/jean-philippe.villemin/extern/hub/hg38/CHIPSEQ_MEAN_EMT/

mv ${OUT}${name}.mean.normalised.bw  ${ARCHIVE_STORE}${name}.mean.normalised.bw

ln -s ${ARCHIVE_STORE}${name}.mean.normalised.bw ${EXTERN}${name}.mean.normalised.bw

#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.rep1.new.bw Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.rep2.new.bw Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.control.new.bw Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.mean.normalised.wig Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.mean.normalised.bw Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.mean.normalised.gz Exists.
#/home/jean-philippe.villemin/data/data/bigWigDeeptools/K27AC_UNT.mean.normalised.png Exists.

exit 0   

