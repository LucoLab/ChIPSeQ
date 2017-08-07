# ChIPSeQ

---

Here we describe the different steps involved in the analysis of ChipSEQ in the laboratory.  

This tutorial goal is to make people understand how the raw data has been process in the team (tools and workflow) and to help them to reproduce analysis.


---

1. Set up Environment
2. Set up directories structure and config file
3. Launch Pipeline


## Set up Environment

---

First you need to check that the following tools are installed on server/computer.

_**Deeptools**_ : Tool to handle deepsequencing files [here](https://deeptools.github.io/)

_**Macs2**_ : Peak Caller [here](https://github.com/taoliu/MACS)

_**wigToBigWig**_ : Include in KentTools suite [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

It's advised to install Conda.   
It will make things easier then for installing tools. (snakemake for example)

_**Conda**_ : [here](https://www.continuum.io/downloads)

_**snakemake**_ : Evil framework for pipeline in python[here](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

_**WiggleTools**_ : _Optionnal_ , can be useful when you want to merge bigwig files and apply some basic numerical operations [here](https://github.com/Ensembl/WiggleTools)

_**SPP**_ :  [here](https://github.com/hms-dbmi/spp)

_**MaSC**_ : [here](http://www.perkinslab.ca/pubs/RPPP2012.html)


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

## Launch Pipeline

In the directory of you condition of interest you must then run : ( Use first --dryrun, -n to see if there is no error in output)

Dry Run : 

```shell
	snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k -n --verbose -p
```

Real Run : 

```shell
	nohup snakemake -s pipeline_chip-seq.Snakefile.py -j 16 -k --verbose -p &
```

Control with htop -u username and look into nohup.out file to see if job is finished.



---



