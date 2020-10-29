# ChIPSeQ

---

Update [here](https://github.com/LucoLab/ChIPSeQ/blob/master/README_update.md)


---

1. Set up Environment - _You don't care !_
2. Set up directories structure and config file
3. Launch Main pipeline
4. Multiqc (marking duplicates and plot number of mapped reads) - _You don't care !_
7. ChipQC
8. DiffBinding and Peak annotation



## Set up directories structure and config file

---

Snakemake pipeline provided by Raoul Raffel @t IGH,Montpellier.

Directory aborescence need to be strictly the same.

One directory per condition. Several histone marks are inside these directories. A _raw_data_ directory where you put fastq files.

Inside this directory, you will set up directories for input raw data and for each replicate of histone marks as followed : 
![ChipSeqDir](https://github.com/ZheFrenchKitchen/pics/blob/master/chipSeqDir.png)

You also have in this directory the script and its config file.

Raw fastq must be formated as follows _file.fq.gz_. 

**NB** : Due to a bug in CTCF and K9AC files. fastq_illumina_filter doesn't work.See ![here](https://github.com/lh3/seqtk/issues/3)
To bypass that, rename the file as _file.filtered.fq.gz_.

**NB1** : If you received different files for one replicate, you need to merge them into one, with a command line like : 

```shell
	cat file1.gz file2.gz > Luco22.fq.gz
```

Then, you need to configure your yaml configuration file. See provided example.

You need to change filepaths,filenames accordingly to your data & environment.


## Launch Main Pipeline

The main pipeline do all alignments for fastq files. It produces also fastqc, wig files for visualisation and finally it will gives us the peaks calling files. Other stuffs are done but you don't care.

Dry Run : 

```shell
	snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k -n --verbose -p
```

Real Run (nohup make the script runs in background. Note the ampersand at the end.): 

```shell
	nohup snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k --verbose -p &
```

Control with htop -u username and look into nohup.out file to see if job is finished.

**NB4** : You can relaunch script if you find something went wrong.It will execute only part that failed. Fastqc will be done again (whereas fastqc analysis has already be done with sucess) because of a bug I didn't manage to correct.



---


<!-- 
## Multiqc

Using all fastqc, you can produce a multiqc wich is a resume of all fastqc produced. Before doint thant you can do a few things to get informations about duplicates.  
You can mark duplicates inside all the bams generated using picard tools.  This is necesseray for the next part when you want to plot a mean signal between replicates whitout replicates reads. 


These commands are optionnal but if you chose to go further, you can then retrieve interesting values that let you plot bargraph with the number of reads mapped, number of reads withtout replicates, etc...see graph below.

```shell
# Flagstat All
find /pathTo/condition -name *.sorted.marked.bam   | xargs -I{} bash -c 'samtools flagstat $0 > ${0/.sorted.marked.bam/.sorted.marked.bash2.txt}' {} \;

# Flagstat concat in one
find /pathTo/condition -name *.sorted.marked.bash2.txt -print -exec cat {} \; > /pathTo/condition/all.flagstat.bash2.txt

# REMOVE DUPLICATES
find /pathTo/condition -name *.sorted.bam   | xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates I=$0 O=${0/.sorted.bam/.sorted.removed.marked.bam} M=${0/.bam/.sorted.removed.marked_dup_metrics.txt} REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT' {} \;

# FLAGSTAT DUPS
find /pathTo/condition -name *.sorted.removed.marked.bam   | xargs -I{} bash -c 'samtools flagstat $0 > ${0/.bam/.sorted.removed.marked.flagstat.bash2.txt}' {}

# Concat in one all flagstat with name of file at the beginning of the block writed
find /pathTo/condition -name *.sorted.removed.marked.flagstat.bash2.txt -print -exec cat {} \; > /pathTo/condition/all.sorted.removed.marked.flagstat.bash2.txt

# Move all fastqc.zip and html for use with multiqc
find pathTo -name *.html | xargs -I{} bash -c 'mv $0 pathTo/FASTQC/$(basename $0)' {} \;
find pathTo -name *.zip | xargs -I{} bash -c 'mv $0 pathTo/FASTQC/$(basename $0)' {} \;
cd pathTo/FASTQC/
multiqc .

```
![Quality](https://github.com/ZheFrenchKitchen/pics/blob/master/K4ME1.curve.png)

To plot this graph, extract interesting values from flagstat ouputs and create the following matrice and save in csv file :

| Sample |	TOTAL |	DUP |	MAP |	TOTAL_WITHOUT_DUP |	MAP_WITHOUT_DUP |
| ---   | --- | --- | --- | --- | --- |
| INPUT_UNT | 77850539 | 12625556 | 75224143 | 65224983 | 62598587 | 
| UNT_2 | 67199665 | 49586720 | 64413504 | 17612945 | 14826784 | 
| T1_1 | 85399142 | 35307778 | 83406602 | 50091364 | 48098824 | 

Then run : 

```shell
Rscript multibarplot.R --file=PathToTheMatriceFile**
```

You can also use this script to plot Number of Raw ChipSeq Peaks, Number of Differentially Bound Histones. (Modify the input file accordingly)
-->
## Normalise tracks

**This step is mandatory if you want to create next bigwig with mean signal without duplicates.**

```shell
#MarkDuplicate
find /pathTo/condition -name *.sorted.bam  | xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates I=$0 O=${0/.sorted.bam/.sorted.marked.bam} M=${0/.bam/.marked_dup_metrics.txt} VALIDATION_STRINGENCY=LENIENT' {} \;
```

I uploaded this part just for record of the methodology I used for the mean track in ucsc browser. It's not automated and path need to be changes, reference files need to be provided, etc, etc....

This has been added to have only one track per mark. 
We used to add each replicates and input in browser. 
Here we are doing one track using mean of replicates removing the input signal signal to normalise.

The script is using a json configuration file. It's also plotting profile around TSS for differentialy expressed gene. So you need to provide bed files for TSS of upregulated genes and TSS of downregulated genes.

**Here we are just interested in creating the mean track, you so can provide fake bed files.**   

```shell
python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/deeptools_countsignal.py --config=/home/jean-philippe.villemin/code/configs/HEATMAP_HISTONE/tss_by_DGE_NODUP/K27AC.json 
```

Other example :  

```shell
python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/deeptools_countsignal.py --config=/home/jean-philippe.villemin/code/configs/HEATMAP_HISTONE/tss_by_DGE_NODUP/POL2_T7.json 
```

_**NB**_ : Lot of things in json you don't care.  
Only things you need to change are the path to bam files.( "bam" & bam-DupRemoved" for Rep1,Rep2,Control)



## ChipQC

First you need to create a samplesheet containing all path to your samples. see( ChipQC ,DiffBind,Doc )

![Quality](https://github.com/ZheFrenchKitchen/pics/blob/master/Listing.png)

You will use this file to create an R object saved on disk that will be used by the next step.
This step is creating a first report with quality controls.


You can use *consensus* methodology where ChIPQC is done using consensus length around peak summit. With this approach, you will also have input signal plotted for interval you defined around peak summits with option -s .  
Option -c is used to say you are using consenus methodology.  
Option -n is used to say what is the name of the Rboject saved.  

```shell
Rscript chipQC.R --file /pathTO/samplesChipMarks.csv -n ALL_MarkedDup -s 500 -c
```

You can find it here : /home/jean-philippe.villemin/code/CHIP-SEQ/R/  

_**NB**_ : I loaded the library at the beginning. It call hg38. (if you work on mouse, you will have to change hard coded stuffs)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

_**NB**_ : Index .bai must be created before. Very important.


![Quality](https://github.com/ZheFrenchKitchen/pics/blob/master/chipQC.png)

## DiffBinding and Peak annotation

Here you will test a variation of binding between condition using R object output (name ALL_MarkedDup in this example) from the step before.  
Option -n is used to say the number of conditions.  
Option -m is used to say which marks from your sample you want to analyze.
Option -p is used to say percentage at least used to retrieve peaks in all you samples. minOverlap (peakset parameter in diffbind) only include peaks in at least this many peaksets in the main binding matrix.If minOverlap is between zero and one, (Don't change 0.99)
peak will be included from at least this proportion of peaksets.  
You will do that for each mark you want to analyse in your sampleshit.

```shell
Rscript chipDiffBind.R --file=/pathTO/ALL_MarkedDup -n 3 -m K27AC -p 0.99
```
**NB** : Don't use the .Rdata extension

Description of output : TODO

---

## Set up Environment

---

First you need to check that the following tools are installed on server/computer.

_**Deeptools**_ : Tool to handle deepsequencing files [here](https://deeptools.github.io/)

_**Macs2**_ : Peak Caller [here](https://github.com/taoliu/MACS)

_**wigToBigWig**_ : Include in KentTools suite [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

It's advised to install Conda.   
It will make things easier then for installing tools. (snakemake for example)

_**Conda**_ : [here](https://www.continuum.io/downloads)

_**snakemake**_ : Evil framework for pipeline in python [here](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

_**WiggleTools**_ : _Optionnal_ , can be useful when you want to merge bigwig files and apply some basic numerical operations [here](https://github.com/Ensembl/WiggleTools)

_**SPP**_ : Need to be there but you don't care [here](https://github.com/hms-dbmi/spp)

_**MaSC**_ : Need to be there but you don't care [here](http://www.perkinslab.ca/pubs/RPPP2012.html)

_**R3.3.1**_ : A lot of packages are used ... It's adviced to read their docs.

_**ChipQC**_ : Control Quality for peaks. Object created will be used then with ChIPQC for Differential binding [here](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html)

_**DiffBind**_ : Differential binding analysis. Which peaks are moving in my samples ? [here](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)

_**ChIPpeakAnno**_ : Anotation of the peaks in a genomic context [here](https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)

_**ChIPseeker**_ : Used for a derived upset peaks graph counting how much peaks are annotated in one or other features. [here](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)

_**rtracklayer**_ : Very usefull to export/import bed to Grange objects. [here](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

_**ComplexHeatmap**_ : Plot heatmap. [here](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

Other packages are mandatory. biomaRt, GenomicFeatures,stringr, reshape , ggplot2, optparse , knitr ... Look inside R scripts to know what to install.

