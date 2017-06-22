#' ---
#' title: "Extract Raw Peaks to bed directly from macs.xls"
#' author: "Jean-Philippe Villemin"
#' date: "June 18, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#'     sample:
#'        value: x
#'     output:
#'        value: x
#' ---

#' ## Description
#' Extract Raw Peaks to bed directly from macs.xls
#' 
#' ## Usage
#' **Rscript extractPeak.R --file=PathToPeaks -s L13S13 -n K27AC_T1**
#'
#' ## Notes
#' None.
#'

#+ setup, include=FALSE
suppressMessages(library(knitr))
suppressMessages(library(ChIPQC))
suppressMessages(library(DiffBind))
suppressMessages(library(rtracklayer))
suppressMessages(library(optparse))

opts_chunk$set(fig.path = 'figure/silk-')

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute File Input Path", metavar="character"),
  make_option(c("-n", "--output"), type="character", default=NULL, help="Name Output", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample name", metavar="character")
); 

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         FUNCTIONS                ###########
######################################################

#' Function to convert MACS output table into a GRanges object :
macs2GRanges <-function(peaks) {
  # generate GRanges object
  myrange <- GRanges(
    seqnames=peaks$chr,
    range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
    strand="*",
    score=0,
    #count=peaks$pileup,
    #score=peaks$X.log10.pvalue.,
    #FE=peaks$fold_enrichment,
    #fdr=peaks$X.log10.qvalue.,
    maxpos=peaks$start+floor(peaks$end-peaks$start)
  )
  return(myrange)
}

######################################################
#########         MAIN                     ###########
######################################################
length(opt)
#' Trick to handle if script is launched directly or by generateDoc
#+ opt-param, include=FALSE
if (length(opt) != 4 ) { 
  opt$file <- params$file
  opt$sample <- params$sample
  opt$output <- params$output
}

filename   = basename(opt$file)

#' Print params
print(opt$file)
print(opt$sample)
print(opt$output)
print(filename)

final_name=paste0(c(opt$sample,opt$output),collapse=".")

chroms=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrY","chrX")

#' Load MACS v2 peaks and convert to GRanges :
macspeaks <- read.table(opt$file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
macspeaks <- macs2GRanges(macspeaks)

#' Clean peaks because browser ucsc doesnt support weird things :
macspeaks_chromOnly <- macspeaks[seqnames(macspeaks) %in% chroms]
export.bed(macspeaks_chromOnly,paste0(c(final_name,"bed"),collapse="."))

