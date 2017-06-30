#!/usr/bin/env bash

echo 'Unmarked standard'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChip.csv -n ALL_unMarkedDup
echo 'Markdup standard'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarks.csv -n ALL_MarkedDup1
echo 'Markdup consenus'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarks.csv -n ALL_MarkedDup2 -c
echo 'Removed Markdup'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipDupRemoved.csv -n ALL_removedDup1
echo 'Removed Markdup Consensus'
#Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipDupRemoved.csv -n ALL_removedDup2  -c
echo 'Markdup consenusK27ME3'
Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarksK27ME3.csv -n ALL_MarkedDupK27ME3_unt2 -c
echo 'Markdup K27ME3'
Rscript /home/jean-philippe.villemin/code/CHIP-SEQ/R/chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarksK27ME3.csv -n ALL_MarkedDupK27ME3_unt
