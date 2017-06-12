#' ---
#' title: "Compare Useq & Macs peaksets by VennDiagram"
#' author: "Jean-Philippe Villemin"
#' date: "June 12, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#'     file2:
#'        value: x
#'     file3:
#'        value: x
#'     outName:
#'        value: x       
#' ---

#' ## Description  
#' Create VennDiagram plots to compare peaks from two differents callers (Macs2 peaks and Useq Peaks).  
#' 
#' One pdf file is create with 3 graphs. One pdf for each histone mark per time point.
#' 
#' ## Usage  
#+ print-command, eval=FALSE
#' **Rscript VennPeaks.R --file=USEQpeaksOutput.xls --file2=MacsRep1.xls --file3=MacsRep2.xls  -o test** 
#' 
#' 
#' 
#+ print-lib, include=FALSE
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(VennDiagram))

#+ print-params, include=FALSE
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="file input from USEQ caller (merged replicate) ", metavar="character"),
  make_option(c("-f2", "--file2"), type="character", default=NULL, help="file input first replicate from MACS caller", metavar="character"),
  make_option(c("-f3", "--file3"), type="character", default=NULL, help="file input second replicate from MACS caller", metavar="character"),
  make_option(c("-o", "--outName"), type="character", default=NULL, help="output name", metavar="character")
); 
parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options


# Usage with zero options passed as parameters in the command line.
# Use parameters set in YAML config at the begining of the file and pass throught render call.
if (length(opt) < 5) { 
  opt$file <- params$file
  opt$file3 <- params$file3
  opt$file2 <- params$file2
  opt$outName <- params$outName
}

print(opt$file)
print(opt$file2)
print(opt$file3)
print(opt$outName)

#' Function to convert MACS output table into a GRanges object :
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
useqpeaks <- read.table(opt$file,
                        skip=1, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                        colClasses=c("NULL", "character", "numeric", "numeric", rep("NULL",7), "numeric",
                                     "numeric", rep("NULL", 7), "numeric", rep("NULL",4))
)
names(useqpeaks) <- c("chr","start","end","summit.start","summit.end", "score")
useqpeaks <- useq2GRanges(useqpeaks[useqpeaks$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])

#+ print-useq, echo=TRUE, eval=FALSE 
head(useqpeaks)

#' Load MACS peaks (only relevant columns) and convert to GRanges :
#+ rep1-macs2, include=TRUE
macspeaks1 <- read.table(opt$file2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
macspeaks1 <- macs2GRanges(macspeaks1[macspeaks1$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])

#+ Print-useq,echo=TRUE
head(macspeaks1)

#' Load MACS peaks and convert to GRanges :
#' (Exactly the same as before for the other replicate)
#+ rep2-macs2
macspeaks2 <- read.table(opt$file3, header=TRUE, sep="\t", stringsAsFactors=FALSE)
macspeaks2 <- macs2GRanges(macspeaks2[macspeaks2$chr %in%  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),])

# /* 
# calculate overlaps between peak sets 
#countOverlaps != findOverlaps, look the doc carrefully for next usage https://rdrr.io/bioc/IRanges/man/findOverlaps-methods.html
#print(length(useqpeaks))
#print(length(macspeaks1))
#length(countOverlaps(macspeaks1, useqpeaks)) # length of macspeaks1 output ?!!
#length(findOverlaps(macspeaks1, useqpeaks,select="all"))  # > length of macspeaks1 output !! query can overlap 2 times in subject, count overlap in subject ?!
#length(macspeaks1 %over% useqpeaks) # finds the ranges in query that overlap any of the ranges in subject. For Ranges or Views objects, it returns a logical vector of length equal to the number of ranges in query
#length(macspeaks1 %within% useqpeaks) #  finds the ranges in query that overlap any of the ranges in subject. For Ranges or Views objects, it returns a logical vector of length equal to the number of ranges in query
#type is within, the query interval must be wholly contained within the subject interval.
#length(subsetByOverlaps(macspeaks1, useqpeaks)) #  returns the subset of query that has an overlap hit with a range in subject */

#' Compute Intersect :
#' 
#+ Overlap, 
UM1=length(subsetByOverlaps(macspeaks1, useqpeaks))
UM2=length(subsetByOverlaps(macspeaks2, useqpeaks))
UM3=length(subsetByOverlaps(macspeaks1, macspeaks2))

#+ Print-Overlap, echo=FALSE, eval=FALSE
print(length(useqpeaks))
print(length(macspeaks1))
print(UM1)

#' Write to file :  
#+ Write, echo=TRUE, eval=FALSE 
pdf(file=paste0(c(opt$outName, "venn.pdf"),collapse="_"))
grid.draw(draw.pairwise.venn( length(useqpeaks), length(macspeaks1), UM1 , 	alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("USeq","firstRep_MACS"),scaled=TRUE))
grid.newpage();
grid.draw(draw.pairwise.venn( length(useqpeaks), length(macspeaks2), UM2 , 	alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("USeq","secondRep_MACS"),scaled=TRUE))
grid.newpage();
draw.pairwise.venn( length(macspeaks1), length(macspeaks2), UM3 , 	alpha = rep(0.3, 2), fill = c("blue", "red"),ext.line.lty = "dashed" ,lty = "blank",category=c("firstRep_MACS2","secondRep_MACS2"),scaled=TRUE)
