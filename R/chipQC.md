ChIP-SEQ QC Analysis
================
Jean-Philippe Villemin
May 22, 2017

-   Resume
    -   Standard
    -   Alternative use (advised)
-   Usage (advised case)
-   Notes
-   Some links that can help

Resume
------

You can apply two strategies depending on the mark you are looking for :

### Standard

You can go with standard method with no options passed in command line.

### Alternative use (advised)

You can use *consensus* methodology where ChIPQ is done using consensus length around peak summit. With this approach, you will also have input signal plotted for interval you defined around peak summits with option -s .Option -c is used to say you are using consenus methodology.

Usage (advised case)
--------------------

Rscript chip.R --file /home/jean-philippe.villemin/CHIPSEQ\_2017\_1\_ALL/samplesChipMarks.csv -n ALL\_MarkedDup -s 500 -c

Notes
-----

This script has been developped as an add-on to what Raoul did first for ChipSeq Analysis.
Raoul snakemake script is doing the alignment with bwa, and peak calling with macs2, some QC control with deeptools and others with spp and MaSc.pl.

Be carefull. You need to mark duplicates with Picard tools if you want duplicates percentage evaluation in your final report.

Some links that can help
------------------------

-   [ChIPQC bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html)
-   [ChIPQC Function Documentation](https://bioconductor.org/packages/release/bioc/manuals/ChIPQC/man/ChIPQC.pdf)
-   [Example of what the report should look like](https://bioconductor.org/packages/release/bioc/vignettes/ChIPQC/inst/doc/ChIPQCSampleReport.pdf)

* * * * *

Needed to avoid memory leaks...By limiting the number of core used by chipQC , you limit memory usage.

``` {.r}
register(MulticoreParam(workers = 10, jobname=opt$name, stop.on.error = TRUE,progressbar = TRUE,log = TRUE, threshold = "INFO",logdir="."), default = TRUE)
```

<!-- register(SerialParam(), default = TRUE)
bpparam() -->



Load a file with experiment details :

``` {.r}
samples <- read.csv(opt$file)
name_file_to_save <- opt$name
```

Do either consensus or standard ChIPQC :

``` {.r}
if (opt$consensus) {
print("Consenus analysis :")
name_file_to_save <- paste0(opt$name,"Consensus",collapse = "_")
name_file_to_save <- paste0(name_file_to_save,opt$summits,collapse = "_")
experiment = ChIPQC(samples,annotation="hg38",consensus=TRUE,bCounts=TRUE,summits=250,chromosomes=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"),blacklist="/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/hg38.blacklist.bed.gz") 
} else
{
print("Standard analysis :")
experiment = ChIPQC(samples,annotation="hg38",chromosomes=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"),blacklist="/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/hg38.blacklist.bed.gz")
}
```

Save RData object on disk to load later for differential affinity binding.

``` {.r}
dbasavefile <- dba.save(experiment, file=name_file_to_save, dir='.', pre='dba_', ext='RData', bMinimize=FALSE) 
```

Create Report.

``` {.r}
ChIPQCreport(experiment,facet=T,lineBy=c("Replicate"),facetBy=c("Condition","Factor"),reportName=paste0(opt$name,"Marks",collapse = "_"), reportFolder=opt$name,colourBy=c("Replicate"))
```
