#' ---
#' title: "Count how many elements there is per condition and events"
#' author: "Jean-Philippe Villemin"
#' date: "Septembrer 1, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#' ---

#' ## Description
#' Count how many elements there is per condition and events
#' 
#' ## Usage
#' **Rscript countOverlap --file=BigOverlapFile.csv**
#'
#' ## Notes
#' None.
#'

#+ setup, include=FALSE
suppressMessages(library(knitr))
suppressMessages(library(optparse))
library(reshape)
library(plyr)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute File Input Path", metavar="character")
); 

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         MAIN                     ###########
######################################################
length(opt)
#' Trick to handle if script is launched directly or by generateDoc
#+ opt-param, include=FALSE
if (length(opt) != 2 ) { 
  opt$file <- params$file
}

filename   = basename(opt$file)
dir_output = dirname(normalizePath(opt$file))

#' Print params
print(opt$file)
print(dir_output)

mydata <- read.csv(file=opt$file, header=FALSE, sep="\t")
count(mydata[,1])

