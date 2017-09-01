#' ---
#' title: "Plot psi vs fold peak"
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
#' Plot labeled point of psi vs fold of chipseq histone peaks
#' 
#' ## Usage
#' **Rscript scatter.R --file=input.csv**
#'
#' ## Notes
#' None.
#'

#+ setup, include=FALSE
suppressMessages(library(knitr))
suppressMessages(library(optparse))
library(reshape)
library(plyr)
library(ggplot2)
library(tidyr)
library(stringr)

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

mydata <- read.csv(file=opt$file, header=FALSE, sep="\t",stringsAsFactors=FALSE)
mydata1 <- data.frame(do.call('rbind', strsplit(as.character(mydata[,5]),'_',fixed=TRUE)))
mydata2 <- data.frame(do.call('rbind', strsplit(as.character(mydata[,1]),'.',fixed=TRUE)))

total <- cbind( mydata1,  mydata2)

clean <- cbind( mydata,  total)

clean <- clean[clean[,16]=="Overlapping",]

filename = gsub(".csv", "", filename)
final=paste0(c(filename,"curve.png"),collapse=".")
final = paste0(c(dir_output,final),collapse="/")
final

head(clean)
clean[,21] <- as.numeric(levels(clean[,21]))[clean[,21]]

# EVENT clean[,31]
# CONDITION clean[,32]
png(file=final,width = 800)
clean_subset <- subset(clean, abs(clean[,12]) >= 0.3 & abs(clean[,21]) >= 1.5)
p <- ggplot(clean, aes(clean[,12], clean[,21]))+ geom_point(aes(colour = factor(clean[,31]),shape = factor(clean[,32])),size=3) + geom_text(data=clean_subset,aes(clean_subset[,12], clean_subset[,21],label=clean_subset[,20],colour = factor(clean_subset[,31])), vjust=2, hjust=1)     + xlab("dpsi")+ ylab("fold peak histone")+labs(color = "Events",shape="Condition")
p
dev.off()

