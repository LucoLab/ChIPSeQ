
# ChIPSeQ -  Update Christmas 2019

## Go on Archive
cd /home/luco/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_CHIPSEQ_FILES/  
mkdir 2019  
cd 2019
 
## Download 
 wget -r -np -bqc -nd -nc --user *** --password *** ftp://blablablah

## Go back to garfield and copy
cd ~/data/data/PROJECT/
mkdir POL2
cd POL2

## Copy 
cp /home/jean-philippe.villemin/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_CHIPSEQ_FILES/2019/RNFC89.fastq.gz .
cp /home/jean-philippe.villemin/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_CHIPSEQ_FILES/2019/RNFC83.fastq.gz .
cp /home/jean-philippe.villemin/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_CHIPSEQ_FILES/2019/RNFC86.fastq.gz .

## Copy 2 files
config.yml
pipeline_chip-seq.Snakefile.py

## Create one dir per replicate. Inside create directory raw_data. Put your fastq.gz

## Rename .fastq.gz to fq.gz

## Modify config.yml

## Run command in directory where you create the directories for each replicate.

```shell
snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k --verbose -p
```

#------------------------------------------------------------
## What next is not automatic, building the wig normalised
#------------------------------------------------------------

## Remove duplicates

```shell
find /pathTo/condition -name *.sorted.bam  | xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates I=$0 O=${0/.sorted.bam/.sorted.marked.bam} M=${0/.bam/.marked_dup_metrics.txt} VALIDATION_STRINGENCY=LENIENT' {} \;
```

## Then move evertyhing to Archive. (You can remove the fastq in raw_data because they should be already backup in archive if you did a copy to garfield. Be sure of that)

## Create the normalised wig (you need to create manually the json file using a given model)

```shell
python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/deeptools_countsignal.py --config=/home/jean-philippe.villemin/code/configs/HEATMAP_HISTONE/tss_by_DGE_NODUP/POL2_T7.json 
```



