#' ---
#' title: "ChIP-SEQ Differential binding Analysis : occupancy & affinity "
#' author: "Jean-Philippe Villemin"
#' date: "June 13, 2017"
#' output:
#'      github_document:
#'        toc: true
#' params:
#'     file:
#'        value: x
#'     numberContrast:
#'        value: x
#'     analyse:
#'        value: x
#'     percent:
#'        value: x
#' ---

#' ## Description
#' All the features from markdown and markdown supported within .Rmd documents, I was able to
#' get from within R scripts.  Here are some that I tested and use most frequently:  
#' 
#' ## Usage
#' **Rscript chipDiffBind.R --file=/home/jean-philippe.villemin/CHIPSEQ_2017_1_ALL/ChIPQC/dba_ALL_MarkedDup2Consensus250 -n 3 -m K27AC**
#'
#' ## Notes
#' Launch command withtout Rdata extensionfor input.
#'


#+ setup, include=FALSE
suppressMessages(library(knitr))
library(ChIPQC)
library(DiffBind)
library(rtracklayer)
suppressMessages(library(optparse))
library(ChIPpeakAnno)
library(biomaRt)
library(GenomicFeatures)
library(ChIPseeker)
library(stringr)


opts_chunk$set(fig.path = 'figure/silk-')

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute File Input Path", metavar="character"),
  make_option(c("-n", "--numberContrast"), type="integer", default=1, help="Nb contrast to check", metavar="number"),
  make_option(c("-m", "--name"), type="character", default=NULL, help="Name of the mark", metavar="character"),
  make_option(c("-p", "--percent"), type="double", default=0.99, help="% of sample to consider per group defaut:0.9", metavar="number"),
  make_option(c("-a", "--analyse"), type="integer", default=1, help="Analysis TimePoint -1 or Mark -2", metavar="integer")
); 

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         FUNCTIONS                ###########
######################################################


#' Function to convert Bed file into a GRanges object :
bed2GRanges <- function(file){
  #https://davetang.org/muse/2015/02/04/bed-granges/
  df <- read.table(file,header=F,stringsAsFactors=F)
  
  if(length(df) > 6){ df <- df[,-c(7:length(df))] }
  if(length(df)<3){ stop("File has less than 3 columns")}
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){ df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)}
  # generate GRanges object
  if(length(df)==3){ gr <- with(df, GRanges(seqnames=chr, range=IRanges(start+1, end))) } 
  else if (length(df)==4){ gr <- with(df, GRanges(seqnames=chr, range=IRanges(start+1, end), id=id)) }
  else if (length(df)==5){  gr <- with(df, GRanges(seqnames=chr, range=IRanges(start+1, end), id=id, score=score))  } 
  else if (length(df)==6){ gr <- with(df, GRanges(seqnames=chr, range=IRanges(start+1, end), id=id, score=score, strand=strand)) }
  
  return(gr)
}

#' Function to convert Bed file into a GRanges object :
gRanges2bed <- function(rangedDataObject){
df <- data.frame(seqnames=seqnames(rangedDataObject),
                 starts=start(rangedDataObject)-1,
                 ends=end(rangedDataObject),
                 names=c(rep(".", length(rangedDataObject))),
                 scores=c(rep(".", length(rangedDataObject))),
                 strands=strand(rangedDataObject))
return(df)
}


######################################################
#########         MAIN                     ###########
######################################################


length(opt)
#' Trick to handle if script is launched directly or by generateDoc
#+ opt-param, include=FALSE
if (length(opt) != 6 ) { 
  opt$file <- params$file
  opt$numberContrast <- params$numberContrast
  opt$name <- params$name
  opt$analyse <- params$analyse
  opt$percent <- params$percent
  
}

filename   = basename(opt$file)
dir_output = dirname(normalizePath(paste0(c(opt$file,"RData"),collapse=".")))

#' Print params
print(opt$file)
print(opt$numberContrast)
print(opt$name)
print(opt$analyse)
print(filename)
print(dir_output)

dir_output = paste0(c(dir_output,opt$name),collapse="/")
print(dir_output)

createDir <- ifelse(!dir.exists(dir_output), dir.create(file.path(dir_output),recursive = TRUE),FALSE)
if(!createDir){print ("Dir already exist")}

dba_chip <- NULL
#Define pdf where graphics will be plotted
filename=paste0(c(filename,opt$name),collapse='.')
final=paste0(c(dir_output,filename),collapse="/")
pdf(file=paste(c(final,"pdf"),collapse="."))
print(final)

filename_analysed = paste0(c(dir_output,"/",filename,"_analysed.RData"),collapse="")
# Name output RData file analysed
print(filename_analysed)

# If analysed RData file doesn"t already exist
if(!file.exists(filename_analysed)){
print("Launch Analysis From Scratch using ChipQC object")
print("LoadChipQCobject")
#' Load dba object : 

dba_chip_init <- dba.load(file=basename(opt$file),dir= dirname(normalizePath(paste0(c(opt$file,"RData"),collapse="."))), pre='', ext='RData')

dba.show(dba_chip_init)

#+ When-you-want-to-recreate-ChipQCreport, echo=FALSE, eval=FALSE
#ChIPQCreport(dba_chip_init,facet=T,lineBy=c("Replicate"),facetBy=c("Condition","Factor"),reportName=paste0(filename,"Marks",collapse = "_"), reportFolder=dir_output,colourBy=c("Replicate"))
#stop()

print("HeatmapOccupancy")
#Correlation heatmap, using occupancy (peak caller score) data
#+ setup, include=TRUE
dba.plotHeatmap(dba_chip_init)

# /* 
# dba.mask(DBA, attribute, value, combine='or', mask, merge='or', bApply=FALSE,
#          peakset, minValue=-1)
# */
# Mask only what you want to check

dba_chip <- dba(dba_chip_init, mask=dba_chip_init@DBA$masks[[opt$name]])
dba_chip$config$fragmentSize <- dba_chip$config$fragmentSize[dba_chip_init@DBA$masks[[opt$name]]]
dba_chip$config$fragmentSize

print("Mask")
dba.show(dba_chip)
names(dba_chip$masks)

dba.plotHeatmap(dba_chip)


type=""
if (opt$analyse==1) {type=DBA_CONDITION
} else { type=DBA_FACTOR}

print("Scatter Overlap")
olap.rate <- dba.overlap(dba_chip,mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

type=""
if (opt$analyse==1) {type=DBA_CONDITION
} else { type=DBA_FACTOR}


# minOverlap only include peaks in at least this many peaksets in the main binding matrix
dba_chip_consensus <- dba.peakset(dba_chip, consensus=c(type),minOverlap=opt$percent)
dba.show(dba_chip_consensus)

dba.plotVenn(dba_chip_consensus,dba_chip_consensus$masks$Consensus)
#dev.off()
consensus_peaks <- dba.peakset(dba_chip_consensus, bRetrieve=TRUE)
# Check bug https://support.bioconductor.org/p/97045/#97051

print("Count")
# /*
#dba.count(DBA, peaks, minOverlap=2, score=DBA_SCORE_TMM_MINUS_FULL, bLog=FALSE,
#          fragmentSize=DBA$config$fragmentSize,
#          summits, filter=0, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
#          mapQCth=DBA$config$mapQCth,
#          filterFun=max,
#          bCorPlot=DBA$config$bCorPlot,
#         bUseSummarizeOverlaps=FALSE, readFormat=DBA_READS_DEFAULT,
#         bParallel=DBA$config$RunParallel)
# minOverlap only include peaks in at least this many peaksets in the main binding matrix
# */
# it first identifies all unique peaks amongst the relevant peaksets and  it merges overlapping peaks
#dba_chip <- dba.count(dba_chip,minOverlap=2)
dba_chip <- dba.count(dba_chip, peaks=consensus_peaks)
#dbasavefile <- dba.save(dba_chip, file=paste0(c(filename,"count"),collapse="_"), dir=dir_output, pre='', ext='RData', bMinimize=FALSE) 
print("Heatmap Consensus Peakset")
dba.show(dba_chip)
dba.plotHeatmap(dba_chip)

print("Contrast")
# /*
# dba.contrast(DBA, group1, group2=!group1, name1="group1", name2="group2",
#              minMembers=3, block, bNot=FALSE,
#              categories=c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT))
# Block could be a useful parameter.. 
# minMembers number of samples at least in one condition
# */

if (opt$analyse==1) {
  dba_chip <- dba.contrast(dba_chip,categories=DBA_CONDITION,minMembers=2)
}  else { 
   dba_chip <- dba.contrast(dba_chip,categories=DBA_FACTOR,minMembers=2)
}

dba.show(dba_chip, bContrast=T)
#dba.plotVenn(dba_chip,contrast=1:opt$numberContrast)

print("Analyse")
# bFullLibrarySize=T method=DBA_DESEQ2
dba_chip <- dba.analyze(dba_chip,bParallel=F) #,bSubControl=FALSE

print("Tcheck if saved")
dbasavefile <- dba.save(dba_chip, file=paste0(c(filename,"analysed"),collapse="_"), dir=dir_output, pre='', ext='RData', bMinimize=FALSE) 

} else {
  
  print("Load Existing DiffBind Analysed Rdata")
  dba_chip <-  dba.load(file=basename(gsub(".RData", "", filename_analysed)), dir=dirname(filename_analysed),pre='', ext='RData') 
  }

print("report")
dba.show(dba_chip)

#The bCalled option adds two columns to the report (Called1 and Called2), one for each group, giving the number of
#samples within the group in which the site was identified as a peak in the original peaksets generated by the peak caller.
#A value of 1 will include all .FDR filter.
# /* dba.report(DBA, contrast, method=DBA$config$AnalysisMethod,
#           th=DBA$config$th, bUsePval=DBA$config$bUsePval,
#           fold=0, bNormalized=TRUE, bFlip=FALSE,
#           bCalled=FALSE, bCounts=FALSE, bCalledDetail=FALSE,
#           bDB, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
#           file,initString=DBA$config$reportInit,ext='csv',DataType=DBA$config$DataType) 
# */
#Retrieve DB sites with FDR < 0.05
#bcalled=True and how many times each site was identified as a peak . This will add a column for
#each group, each indicating the number of samples in the group identified as a peak in the original peaksets.
#bcounts=True normalized counts
#bCalledDetail peak caller status should be included for each sample (if available). Columns are added for each sample in the first group, followed by
#columns for each sample in the second group.
#bFlip

#dba_chip <-  dba.load(file=paste0(c(filename,"K27AC_count_with_analysis"),collapse="."), dir=paste0(c(dir_output,"/K27AC"),collapse=""), pre='', ext='RData') 

txdb <- loadDb("/home/jean-philippe.villemin/gencode.v25.annotation.sqlite")
edb = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="jul2016.archive.ensembl.org")

# Test each contrast
for (contrast in 1:opt$numberContrast) {

  print(dba_chip$contrast[[contrast]]$name1)
  print(dba_chip$contrast[[contrast]]$name2)
  
  contrast_name=paste0(c(dba_chip$contrast[[contrast]]$name1,"_vs_",dba_chip$contrast[[contrast]]$name2),collapse="")
  
  filename2=paste0(c(filename,contrast_name),collapse=".")
  print(filename2)
  # Export classical csv
  dba_chip.DB <-dba.report(dba_chip,contrast=contrast,bCalled=TRUE,bCounts=TRUE,bCalledDetail=TRUE,ext='tsv',file=filename2,initString="")
  
  
  if(is.null(dba_chip.DB)){next}
  print(paste0(c(getwd(),"/","_",filename2,".tsv"),collapse=""))
  print(paste0(c(dirname(filename_analysed),'/',filename2,".tsv"),collapse=""))
  # Mv file from current directory of the script to output
  
  file.copy(from = paste0(c(getwd(),"/","_",filename2,".tsv"),collapse=""), to = paste0(c(dirname(filename_analysed),'/',filename2,".tsv"),collapse=""))
  
  file.remove(paste0(c(getwd(),"/","_",filename2,".tsv"),collapse=""))
  
  chroms=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrY","chrX")

  
  # Grep DBA_DATA_GRANGES
  dba_chip_Grange <-dba.report(dba_chip,contrast=contrast,DataType=DBA_DATA_GRANGES)

  export.bed(dba_chip_Grange[seqnames(dba_chip_Grange) %in% chroms], paste0(c(final,contrast_name,"diffPeaks.bed"),collapse="."))
  
  print("Annote peaks")
  # Annote peaks
  dba_chip_Grange.anno <- annotatePeakInBatch(dba_chip_Grange, 
                                              featureType="TSS", 
                                              AnnotationData =  toGRanges(txdb), 
                                              output = "both",  # will output all the nearest features, in addition, will output any features that overlap the peak that is not the nearest features
                                              select="all", # may return multiple overlapping peaks
                                              maxgap = 0L, 
                                              ignore.strand = FALSE, 
                                              PeakLocForDistance = "middle", # Specify the location of peak for calculating distance,i.e., middle means using middle of the peak to calculate distance to feature
                                              FeatureLocForDistance = "TSS", # TSS means using start of feature when feature is on plus strand and using end of feature when feature is on minus strand
                                              bindingRegion = NULL) # Specifying the criteria to associate peaks with annotation .see annoPeaks. Use with FeatureLocForDistance . Compute region to look for around PeakForDistance. Argument Can be c(-5000,+5000)
  #From doc possible values  :
  #insideFeature: upstream peak resides upstream of the feature; downstream: peak resides downstream of the feature; inside: peak resides inside the feature; overlapStart: peak overlaps with the start of the feature; overlapEnd: peak overlaps with the end of the feature; includeFeature: peak include the feature entirely
  #distancetoFeature: distance to the nearest feature such as transcription start site. By default, the distance is calculated as the distance between the start of the binding site and the TSS that is the gene start for genes located on the forward strand and the gene end for genes located on the reverse stran
  #fromOverlappingOrNearest: nearest: indicates this feature's start (feature's end for features at minus strand) is closest to the peak start;
  #shortestDistance: The shortest distance from either end of peak to either end the feature.
  
  # Remove version of ensembl gene ENSXXXXX.X to use with addGeneIDs
  dba_chip_Grange.anno$feature <- str_match(dba_chip_Grange.anno$feature,"^(\\w+)\\.([0-9]+)")[,2]
  
  # Add Symbol
  dba_chip_Grange.fullanno <- addGeneIDs(dba_chip_Grange.anno, mart=edb, feature_id_type="ensembl_gene_id", IDs2Add= c("hgnc_symbol","gene_biotype"), silence=FALSE)
  # Sometime feature is NA. WTF. Remove these lines.
  dba_chip_Grange.fullanno <- dba_chip_Grange.fullanno[which(!is.na(dba_chip_Grange.fullanno$feature))]
  # Transform pvalue using log10
  dba_chip_Grange.fullanno$`p-value` <- -log10(dba_chip_Grange.fullanno$`p-value`)
  # If Symbol not found and value empty, set either Ensembl ID 
  dba_chip_Grange.fullanno$hgnc_symbol[dba_chip_Grange.fullanno$hgnc_symbol == ""] <- dba_chip_Grange.fullanno$feature
  # remove _ for protein_coding to avoid splitting problem next
  dba_chip_Grange.fullanno$gene_biotype <- gsub("_", "", dba_chip_Grange.fullanno$gene_biotype)
  
  # Create ID that will be exported in bed
  names(dba_chip_Grange.fullanno) <- paste(dba_chip_Grange.fullanno$peak,dba_chip_Grange.fullanno$fromOverlappingOrNearest,dba_chip_Grange.fullanno$insideFeature,dba_chip_Grange.fullanno$feature,dba_chip_Grange.fullanno$gene_biotype,dba_chip_Grange.fullanno$hgnc_symbol,dba_chip_Grange.fullanno$Fold,dba_chip_Grange.fullanno$`p-value`,dba_chip_Grange.fullanno$feature_strand,dba_chip_Grange.fullanno$distancetoFeature,dba_chip_Grange.fullanno$shortestDistance, sep="_")
  
  head(dba_chip_Grange.fullanno)
  
  # Filter out stange seqnames and export
  clean <- dba_chip_Grange.fullanno[seqnames(dba_chip_Grange.fullanno) %in% chroms]
  # Export bed
  export.bed(clean, paste0(c(final,contrast_name,"bed"),collapse="."))

  #look at the binding affinity differences
  dba.plotMA(dba_chip,contrast=contrast)
  if (opt$analyse==1){ dba.plotPCA(dba_chip, DBA_CONDITION, label=DBA_CONDITION,contrast=contrast) } 
  else { dba.plotPCA(dba_chip, DBA_FACTOR, label=DBA_FACTOR,contrast=contrast) } 
  dba.plotPCA(dba_chip, attributes=c(DBA_FACTOR,DBA_CONDITION),label=DBA_ID,contrast=contrast)
  
  #Correlation heatmap, using only significantly differentially bound sites.
  dba.plotHeatmap(dba_chip, contrast=contrast,correlations=TRUE)
  dba.plotHeatmap(dba_chip, contrast=contrast, correlations=FALSE)
  dba.plotHeatmap(dba_chip,contrast=contrast, correlations=FALSE, scale="row")
  dba.plotBox(dba_chip,contrast=contrast)


  # Annote around TSS with ChipSeqSeeker
  #peakAnno <- annotatePeak(clean, tssRegion=c(-3000, 3000), TxDb=txdb)
  
  #png(file=paste0(c(final,contrast,"camenbert","png"),collapse="."))
  #plotAnnoPie(peakAnno)
  #dev.off()
  
  #png(file=paste0(c(final,contrast,"upset","png"),collapse="."))
  #upsetplot(peakAnno)
  #dev.off()
  
  
}
dev.off()


#/* LonG LonG LonG Description of what report can do
#if neither bDB or bNotDB is set to TRUE, a report dataframe, GRanges, or RangedData object is
#returned, with a row for each binding site within the thresholding parameters, and the following
#columns:
#  Chr Chromosome of binding site
#Start Starting base position of binding site
#End End base position of binding site
#Conc Concentration – mean (log) reads across all samples in both groups
#Conc_group1 Group 1 Concentration – mean (log) reads across all samples first group
#Conc_group2 Group 2 Concentration – mean (log) reads across all samples in second group
#Fold Fold difference – mean fold difference of binding affinity of group 1 over group
#2 (Conc1 - Conc2). Absolute value indicates magnitude of the difference, and
#sign indicates which one is bound with higher affinity, with a positive value
#indicating higher affinity in the first group
#p-value p-value calculation – statistic indicating significance of difference (likelihood
#difference is not attributable to chance)
#FDR adjusted p-value calculation – p-value subjected to multiple-testing correction
#If bCalled is TRUE and caller status is available, two more columns will follow:
#  Called1 Number of samples in group 1 that identified this binding site as a peak
#Called2 Number of samples in group 2 that identified this binding site as a peak
#If bCounts is TRUE, a column will be present for each sample in group 1, followed by each sample
#in group 2. The Sample ID will be used as the column header. This column contains the read counts
#for the sample.
#If bCalledDetail is TRUE, a column will be present for each sample in group 1, followed by each
#sample in group 2. The Sample ID will be used as the column header. This column contains a "+"
#to indicate for which sites the sample was called as a peak, and a "-" if it was not so identified.
#If bDB or bNotDB is set to TRUE, a special DBA object is returned,containing peaksets based on sites
#determined to be differentially bound (or not) as specified using the bDB, bNotDB, bGain, bLoss,
#and bAll parameters. In this DBA object, the Tissue value will specify the direction of the change
#(Gain for positive fold changes, Loss for negative fold changes, and All for any fold change).
#The Factor value specifies if the peaks are differentially bound (DB) or not (!DB). The Condition
#value specifies the analysis method (e.g. edgeR)
#*/
