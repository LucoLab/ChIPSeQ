ChIP-SEQ QC Analysis
================
Jean-Philippe Villemin
May 22, 2017

ChipQC <br> \* h1 \* h2 \* h3 <br> \# Control Quality for ChipSeq You can apply two strategies depending on the mark you are looking for :
First <br> \#\# Notes : This script has been developped as an add-on to what Raoul did first for ChipSeq Analysis. ? <!-- Todo find Raoul surname... -->
Raoul script are doing the alignment with bwa, and variant calling with macs2 , some QC control with deeptools and other with spp and MaSc.pl. <br> Be carefull. You need to mark duplicates with Picard tools if you want duplicates percentage evaluation in your final report.
\#\#\# Some links that can help : \#\#\# \* [ChIPQC bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html) \* [ChIPQC Function Documentation](https://bioconductor.org/packages/release/bioc/manuals/ChIPQC/man/ChIPQC.pdf) \* [Example of what the report should look like](https://bioconductor.org/packages/release/bioc/vignettes/ChIPQC/inst/doc/ChIPQCSampleReport.pdf)

``` {.r}
print(opt$file)
```

    ## [1] "/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChip.csv"

``` {.r}
print(opt$name)
```

    ## [1] "Test"

``` {.r}
#names(registered())
#print('Registred')
#registered()
```

Needed to avoid memory leaks...By limiting the number of core used by chipQC , you limit memory usage.

``` {.r}
register(MulticoreParam(workers = 16, jobname=opt$name, stop.on.error = TRUE,progressbar = TRUE,log = TRUE, threshold = "INFO",logdir="."), default = TRUE)
#register(SerialParam(), default = TRUE)
#bpparam()
```

Load file with experiment details

``` {.r}
samples <- read.csv(opt$file)
```
