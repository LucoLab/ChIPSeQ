#' ---
#' title: "Extract peakset from USEQ to write peaks"
#' author: "Jean-Philippe Villemin"
#' date: "June 12, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#'     outNameBed:
#'        value: x
#' ---

#' ## Description  
#' Output to bed file peaksets from USeq PeakCaller
#' 
#' ## Usage  
#' 
#+ print-command, eval=FALSE
#' **Rscript USEQ.R --file=UseqPeaksoutput.xls -o test**
#' 
#'  

#+ printlib, include=FALSE
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="file input", metavar="character"),
  make_option(c("-o", "--outNameBed"), type="character", default=NULL, help="Bed output name", metavar="character")
); 
parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

# Usage with zero options passed as parameters in the command line.
# Use parameters set in YAML config at the begining of the file and pass throught render call.

if (length(opt) != 3) { 
opt$file <- params$file
opt$outNameBed <- params$outNameBed
}

print(opt$file)
print(opt$outNameBed)

#' Function to convert USeq output table into a GRanges object :
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

#' Load USeq peaks (only relevant columns) and convert to GRanges :
useqpeaks <- read.table(opt$file, skip=1, header=FALSE, sep="\t", stringsAsFactors=FALSE, colClasses=c("NULL", "character", "numeric", "numeric", rep("NULL",7), "numeric", "numeric", rep("NULL", 7), "numeric", rep("NULL",4)))
names(useqpeaks) <- c("chr","start","end","summit.start","summit.end", "score")
useqpeaks <- useq2GRanges(useqpeaks)

#' Pass to dataFrame and select columns of interest to construct your bed format :
#+ data.frame, include=TRUE
df <- data.frame(seqnames=seqnames(useqpeaks),
                 starts=start(useqpeaks)-1,
                 ends=end(useqpeaks),
                 names=c(rep(".", length(useqpeaks))),
                 scores=elementMetadata(useqpeaks)$score,
                 strands=strand(useqpeaks))

#' Write to file :
write.table(df, file=paste0(c(opt$outNameBed,".bed"),collapse=""), quote=F, sep="\t", row.names=F, col.names=F)