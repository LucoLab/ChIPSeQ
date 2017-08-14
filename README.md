# ChIPSeQ

---

Here we describe the different steps involved in the analysis of ChipSEQ in the laboratory.  

This tutorial goal is to make people understand how the raw data has been processed in the team (tools and workflow) and to help them to reproduce analysis.


---

1. Set up Environment
2. Set up directories structure and config file
3. Launch Main pipeline
4. Multiqc (marking duplicates and plot number of mapped reads)
7. ChipQC
8. DiffBinding and Peak annotation

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

Other packages are mandatory. biomaRt, GenomicFeatures ,stringr, reshape , ggplot2, optparse , knitr ... Look inside R scripts to know what to install.

## Set up directories structure and config file

---

Analyse is based on snakemake pipeline provided by Raoul Raffel @t IGH,Montpellier.

Directory aborescence need to be strictly the same.

One directory per condition. Several histone marks are inside these directories.

Inside this directory, you will set up directories for input raw data and for each replicate of histone marks as followed : 
![ChipSeqDir](https://github.com/ZheFrenchKitchen/pics/blob/master/chipSeqDir.png)

As you can see in the picture, you also have in this directory the script and its config file.

In directories for the replicates only (not input), you need to place your raw data in _raw_data_ directory. Raw data must be formated as follows _file.fq.gz_.

**NB1** : If you received different files for one replicate, you need to merge them into one, with a command line like : 

```shell
	cat file1.gz file2.gz > Luco22.fq.gz
```

Then, you need to configure your yaml configuration file. See provided example.

You need to change filepaths,filenames accordingly to your data & environment.

**NB2** : In script, there is hard coded path, you need to manually change in code "--tempdir /home/jean-philippe.villemin/tmp_deeptools/ ". Create a dir wherever you want and change in the code.
This is used by deeptools to write temporary files...

**NB3** : You need to download Rscripts in dependencies dir from github. Set up paths in configfile. Same has to be done for the bedfile with exons used to do heatmaps and other QC stuffs. Exons were randomly selected. You can provide your own list.

## Launch Main Pipeline

The main pipeline do all alignments for fastq files. It produces also fastqc, wig files for visualisation and finally it will gives us the peaks calling files. Other stuffs are done but are less important to continue the tutorial.

In the directory of you condition of interest you must then run : ( Use first --dryrun, -n to see if there is no error in output)

Dry Run : 

```shell
	snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k -n --verbose -p
```

Real Run (nohup make the script runs in background.Note the ampersand at the end.): 

```shell
	nohup snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k --verbose -p &
```

Control with htop -u username and look into nohup.out file to see if job is finished.

**NB4** : You can relaunch script if you find something went wrong.It will execute only part that failed. Fastqc will be done again (whereas fastqc analysis has already be done with sucess) because of a bug I didn't manage to correct.

At the end of this part, you should have bam, bigwig and file containing peaks calling.

---

## Multiqc

Using all fastqc, you can produce a multiqc wich is a resume of all fastqc produced. Before doint thant you can do a few things to get informations about duplicates.  
You can mark duplicates inside all the bams generated using picard tools.  
You can retrieve interesting values that let you plot bargraph with the number of reads mapped, number of reads withtout replicates, etc...see graph below.

```shell
#MarkDuplicate
find /pathTo/condition -name *.sorted.bam  | xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates I=$0 O=${0/.sorted.bam/.sorted.marked.bam} M=${0/.bam/.marked_dup_metrics.txt} VALIDATION_STRINGENCY=LENIENT' {} \;

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

## ChipQC

First you need to create a samplesheet containing all path to your samples. see( ChipQC ,DiffBind,Doc )

![Quality](https://github.com/ZheFrenchKitchen/pics/blob/master/Listing.png)

You will use this file to create an R object saved on disk that will be used by the next step.
This step is creating a first report with quality controls.


You can use *consensus* methodology where ChIPQC is done using consensus length around peak summit. With this approach, you will also have input signal plotted for interval you defined around peak summits with option -s .  
Option -c is used to say you are using consenus methodology.  
Option -n is used to say what is the name of the Rboject saved.  

```shell
Rscript chip.R --file /pathTO/samplesChipMarks.csv -n ALL_MarkedDup -s 500 -c
```

## DiffBinding and Peak annotation

Here you will test a variation of binding between condition using output from the step before.  
Option -n is used to say the number of conditions.  
Option -m is used to say which marks from your sample you want to analyze.

```shell
Rscript chipDiffBind.R --file=/pathTO/dba_ALL_MarkedDup2Consensus250 -n 3 -m K27AC
```



