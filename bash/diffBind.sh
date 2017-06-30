#!/usr/bin/env bash

# NB : K27ME3 k27 is done using 0.1 . It means at least one peaks in one replicate or the other is taken into account. 

echo 'K4ME1'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_MarkedDup2Consensus250 -n 3 -m K4ME1

echo 'K27ME3'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_MarkedDup2Consensus250 -n 3 -m K27ME3

echo 'K27AC'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_MarkedDup2Consensus250 -n 3 -m K27AC

echo 'K27AC removedDup'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_removedDup2Consensus250 -n 3 -m K27AC
echo 'K27ME3 removedDup'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_removedDup2Consensus250 -n 3 -m K27ME3
echo 'K4ME1 removedDup'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_removedDup2Consensus250 -n 3 -m K4ME1
