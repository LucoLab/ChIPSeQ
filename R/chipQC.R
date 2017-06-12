#' ---
#' title: "ChIP-SEQ QC Analysis"
#' author: "Jean-Philippe Villemin"
#' date: "June 7, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     name:
#'        value: x
#'     file:
#'        value: x
#'     consensus:
#'        value: x
#'     summits:
#'        value: x       
#' ---


#' ## Resume  
#' You can apply two strategies depending on the mark you are looking for :  
#' 
#' ### Standard  
#' 
#' You can go with standard method with no options passed in command line.  
#' 
#' ### Alternative use (advised)  
#' 
#' You can use *consensus* methodology where ChIPQ is done using consensus length around peak summit. With this approach, you will also have input signal plotted for
#' interval you defined around peak summits with option -s .Option -c is used to say you are using consenus methodology.
#' 
#' ## Usage (advised case)
#+ usage, echo=TRUE, eval=FALSE
#' **Rscript chip.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarks.csv -n ALL_MarkedDup -s 500 -c**
#'  
#' ## Notes  
#' 
#' This script has been developped as an add-on to what Raoul did first for ChipSeq Analysis.  
#' Raoul snakemake script is doing the alignment with bwa, and peak calling with macs2, some QC control with deeptools and others with spp and MaSc.pl.  
#' 
#' Be carefull. You need to mark duplicates with Picard tools if you want duplicates percentage evaluation in your final report.
#' 
#' ## Some links that can help  
#' 
#'    * [ChIPQC bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html)  
#'    * [ChIPQC Function Documentation](https://bioconductor.org/packages/release/bioc/manuals/ChIPQC/man/ChIPQC.pdf)  
#'    * [Example of what the report should look like](https://bioconductor.org/packages/release/bioc/vignettes/ChIPQC/inst/doc/ChIPQCSampleReport.pdf)  
#'     
#' 
#' ---


#+ setup, include=FALSE
suppressMessages(library(knitr))
library(ChIPQC)
library(BiocParallel)
suppressMessages(library(optparse))
opts_chunk$set(fig.path = 'figure/silk-')

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="File input", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, help="Name output ", metavar="character"),
  make_option(c("-c", "--consensus"),action="store_true", default=FALSE,help="Chose if you want to use consensus mode"),
  make_option(c("-s", "--summits"), action="store", default=250,help="Length up and downstream the peak when consensus is used only [default %default] * 2") 
  
); 

# Rscript chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarks.csv -n ALL_MarkedDup -s 500 -c
# Rscript chipQC.R --file /home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/samplesChipMarks.csv -n ALL_MarkedDup

parser = OptionParser(usage = "%prog [options] ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options


print(length(params))
print(length(opt))

# Usage with zero options passed as parameters in the command line.
# Use parameters set in YAML config at the begining of the file and pass throught render call.
if (length(opt) == 3) { 
opt$file <- params$file
opt$name <- params$name
opt$consensus <- params$consensus
opt$summits <- params$summits
}

#+ print-params, include=FALSE
print(opt$file)
print(opt$name)
print(opt$consensus)
print(opt$summits)

# /* names(registered())
# print('Registred')
# registered() */

#' Needed to avoid memory leaks...By limiting the number of core used by chipQC , you limit memory usage.
#+ register-cores,  include=TRUE
register(MulticoreParam(workers = 10, jobname=opt$name, stop.on.error = TRUE,progressbar = FALSE,log = FALSE, threshold = "INFO",logdir="."), default = TRUE)
#' <!-- register(SerialParam(), default = TRUE)
#' bpparam() -->

#+ load-experiment, include=TRUE
#' Load a file with experiment details :
samples <- read.csv(opt$file)
name_file_to_save <- opt$name


#' Do either consensus or standard ChIPQC :
#+ QC-experiment-consensus, echo=TRUE, eval=FALSE
if (opt$consensus) {
print("Consenus analysis :")
name_file_to_save <- paste0(opt$name,"Consensus",collapse = "_")
name_file_to_save <- paste0(name_file_to_save,opt$summits,collapse = "_")
experiment = ChIPQC(samples,annotation="hg38",consensus=TRUE,bCounts=TRUE,summits=250,chromosomes=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),blacklist="/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/hg38.blacklist.bed.gz") 
} else
{
print("Standard analysis :")
experiment = ChIPQC(samples,annotation="hg38",chromosomes=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),blacklist="/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/hg38.blacklist.bed.gz")
}

#' Save RData object on disk to load later for differential affinity binding.
#+ save-object, echo=TRUE, eval=FALSE
dbasavefile <- dba.save(experiment, file=name_file_to_save, dir='.', pre='dba_', ext='RData', bMinimize=FALSE) 

#' Create Report.
#+ do-report, echo=TRUE, eval=FALSE
ChIPQCreport(experiment,facet=T,lineBy=c("Replicate"),facetBy=c("Condition","Factor"),reportName=paste0(opt$name,"Marks",collapse = "_"), reportFolder=opt$name,colourBy=c("Replicate"))
