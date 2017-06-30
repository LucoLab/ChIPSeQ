#!/usr/bin/env bash


#find /home/jean-philippe.villemin/data/data/workspace/beds/exon/met -name *.bed | xargs -I{} bash -c 'Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/localisation.R --file=$0' {} \;

#find /home/jean-philippe.villemin/extern/secure/listing/intersect -name *7.hg38.vs.ExonEMT.bed | xargs -I{} bash -c 'Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/localisation.R --file=$0' {} \;

#find /home/jean-philippe.villemin/RNASEQ_2017_RESULTS/IRFINDER/ -name *hg38.sorted.bed | xargs -I{} bash -c 'Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/localisation.R --file=$0' {} \;


find /home/jean-philippe.villemin/RNASEQ_2017_RESULTS/RMATS/19_5_2017__14_34_10/ -name *.bed | xargs -I{} bash -c 'Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/localisation.R --file=$0' {} \;


find /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/RawPeaks/ -name *.bed | xargs -I{} bash -c 'Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/localisation.R --file=$0' {} \;

