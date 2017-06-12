Compare Useq & Macs peaksets by VennDiagram
================
Jean-Philippe Villemin
June 12, 2017

-   Description
-   Usage

Description
-----------

Create VennDiagram plots to compare peaks from two differents callers (Macs2 peaks and Useq Peaks).

One pdf file is create with 3 graphs. One pdf for each histone mark per time point.

Usage
-----

**Rscript VennPeaks.R --file=USEQpeaksOutput.xls --file2=MacsRep1.xls --file3=MacsRep2.xls -o test**

Function to convert MACS output table into a GRanges object :

``` {.r}
macs2GRanges <-function(peaks) {
  # generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
    strand="*",
    count=peaks$pileup,
    score=peaks$X.log10.pvalue.,
    FE=peaks$fold_enrichment,
    fdr=peaks$X.log10.qvalue.,
    maxpos=peaks$start+floor(peaks$end-peaks$start)
  )
  return(myrange)
}
```

Function to convert USeq output table into a GRanges object :

``` {.r}
useq2GRanges <-function(peaks) {
  # generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
    strand="*",
    score=peaks$score,
    # use midpoint of summit bin as maxpos
    maxpos=peaks$summit.start+floor(peaks$summit.end-peaks$summit.start)
  )
  
  return(myrange)
}
```

Load USeq peaks (only relevant columns) and convert to GRanges :

``` {.r}
useqpeaks <- read.table(opt$file,
                        skip=1, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                        colClasses=c("NULL", "character", "numeric", "numeric", rep("NULL",7), "numeric",
                                     "numeric", rep("NULL", 7), "numeric", rep("NULL",4))
)
names(useqpeaks) <- c("chr","start","end","summit.start","summit.end", "score")
useqpeaks <- useq2GRanges(useqpeaks[useqpeaks$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])
```

``` {.r}
head(useqpeaks)
```

Load MACS peaks (only relevant columns) and convert to GRanges :

``` {.r}
macspeaks1 <- read.table(opt$file2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
macspeaks1 <- macs2GRanges(macspeaks1[macspeaks1$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])
```

``` {.r}
head(macspeaks1)
```

    ## GRanges object with 6 ranges and 5 metadata columns:
    ##              seqnames         ranges strand |     count     score
    ##                 <Rle>      <IRanges>  <Rle> | <numeric> <numeric>
    ##   chr1:14021     chr1 [14021, 14532]      * |     17.37   6.87292
    ##   chr1:22215     chr1 [22215, 22576]      * |     14.77   4.35855
    ##   chr1:25743     chr1 [25743, 26520]      * |     12.34   3.44742
    ##   chr1:27273     chr1 [27273, 28581]      * |     18.57   8.56743
    ##   chr1:29491     chr1 [29491, 30589]      * |      17.6   7.90214
    ##   chr1:31206     chr1 [31206, 31342]      * |     12.38     3.715
    ##                     FE       fdr    maxpos
    ##              <numeric> <numeric> <numeric>
    ##   chr1:14021   3.79911   5.20551     14532
    ##   chr1:22215   2.88476   2.87945     22576
    ##   chr1:25743   2.61386   2.05308     26520
    ##   chr1:27273    4.2871   6.79881     28581
    ##   chr1:29491   4.21992   6.16576     30589
    ##   chr1:31206    2.7951   2.29307     31342
    ##   -------
    ##   seqinfo: 23 sequences from an unspecified genome; no seqlengths

Load MACS peaks and convert to GRanges : (Exactly the same as before for the other replicate)

``` {.r}
macspeaks2 <- read.table(opt$file3, header=TRUE, sep="\t", stringsAsFactors=FALSE)
macspeaks2 <- macs2GRanges(macspeaks2[macspeaks2$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])
```

Compute Intersect :

``` {.r}
UM1=length(subsetByOverlaps(macspeaks1, useqpeaks))
UM2=length(subsetByOverlaps(macspeaks2, useqpeaks))
UM3=length(subsetByOverlaps(macspeaks1, macspeaks2))
```

Write to file :

``` {.r}
pdf(file=paste0(c(opt$outName, "venn.pdf"),collapse="_"))
grid.draw(draw.pairwise.venn( length(useqpeaks), length(macspeaks1), UM1 ,  alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("USeq","firstRep_MACS"),scaled=TRUE))
grid.newpage();
grid.draw(draw.pairwise.venn( length(useqpeaks), length(macspeaks2), UM2 ,  alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("USeq","secondRep_MACS"),scaled=TRUE))
grid.newpage();
draw.pairwise.venn( length(macspeaks1), length(macspeaks2), UM3 ,   alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("firstRep_MACS2","secondRep_MACS2"),scaled=TRUE)
```
