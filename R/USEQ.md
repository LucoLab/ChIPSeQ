USEQ write to file peaksets
================
Jean-Philippe Villemin
June 12, 2017

-   Description
-   Usage

Description
-----------

Output to bed file peaksets from USeq PeakCaller

Usage
-----

**Rscript USEQ.R --file=UseqPeaksoutput.xls -o test**

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
useqpeaks <- read.table(opt$file, skip=1, header=FALSE, sep="\t", stringsAsFactors=FALSE, colClasses=c("NULL", "character", "numeric", "numeric", rep("NULL",7), "numeric", "numeric", rep("NULL", 7), "numeric", rep("NULL",4)))
names(useqpeaks) <- c("chr","start","end","summit.start","summit.end", "score")
useqpeaks <- useq2GRanges(useqpeaks)
```

Pass to dataFrame and select columns of interest to construct your bed format :

``` {.r}
df <- data.frame(seqnames=seqnames(useqpeaks),
                 starts=start(useqpeaks)-1,
                 ends=end(useqpeaks),
                 names=c(rep(".", length(useqpeaks))),
                 scores=elementMetadata(useqpeaks)$score,
                 strands=strand(useqpeaks))
```

Write to file :

``` {.r}
write.table(df, file=paste0(c(opt$outNameBed,".bed"),collapse=""), quote=F, sep="\t", row.names=F, col.names=F)
```
