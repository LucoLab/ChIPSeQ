#' ---
#' title: "Plot Genomic Coverage per % Reads Sampled "
#' author: "Jean-Philippe Villemin"
#' date: "June 29, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#' ---

#' ## Description
#' Extract Raw Peaks to bed directly from macs.xls
#' 
#' ## Usage
#' **Rscript coverageByReadSamp.R --file=PathToPeaks**
#'
#' ## Notes
#' None.
#'

#+ setup, include=FALSE
suppressMessages(library(knitr))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
library(reshape)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute File Input Path", metavar="character")
); 

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         FUNCTIONS                ###########
######################################################



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

mydata <- read.csv(file=opt$file, header=TRUE, sep=",")
#df <- read.table(,header = TRUE,sep = ",")
mydata
dfm <- melt(mydata[,c('Sample','TOTAL','MAP','DUP','TOTAL_WITHOUT_DUP','MAP_WITHOUT_DUP')],id.vars = 1)

filename = gsub(".csv", "", filename)
final=paste0(c(filename,"curve.png"),collapse=".")
final = paste0(c(dir_output,final),collapse="/")
final
png(file=final,width = 800)
#ggplot(mydata, aes(x=mydata$V2,y=mydata$V3,color=mydata$V1))  + labs(x = "%Sampling",y="Genomic Coverage bp",title="Genomic Coverage per % Reads Total Sampled") + geom_point(shape=1) + geom_line()
ggplot(dfm, aes(x = Sample,y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge")
dev.off()
#


