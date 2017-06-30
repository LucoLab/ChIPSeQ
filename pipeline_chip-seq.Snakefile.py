"""
@version: 1.0
@author: Raoul RAFFEL

This pipeline is under developpment but main features to analyse ChIP-seq
is implemented. 

Tools used :
    Quality Control   : 
            - fastQC (sequencing quality)
            - SPP & MaSC (cross strand correlation : fraglength)
            - deeptool : 
                - bamFigerprint : broad/narrow peaks
                - correlation : IP replicate good correlation
                                IP/input bad correlation
                - GC bias (doesn't work)
                TODO : 
                    - depptools PCA

    Mapping           : BWA.

    bam => bigwig     : deeptools bamCoverage & bigwigCompare

    Peaks calling     : MACS2 (broad/narrow peaks)
                            TODO:
                                SPP
                                coda
                                

    intersect peaks   : awk & bedtools & python 

    heatmap & profile : deeptools computeMatrix, plotHeatmap, plotProfile


Deeptools :
     make bigwig files normalised on sequencing depth (RPKM)
     make bigwig normalised on input (ratio)
     make heatmap/profile of Start, End, center and body of region files (annotations in bed format)

use bedtools and awk to intersect peaks and features (direct_target_peaks rule)

use config.yml/(config.json) to choose what kind of heatmap/profile to make (TSS, TES, body, center)

/!\ Note : json format is deprecated. You can still use it but it's more readable to use yaml. YAML allow comment (#) 
"""
import os
import re

import subprocess
# get version of a tools i.e.:
# i.e.: subprocess.getoutput("bamCoverage --version")

# import config file 
if os.path.isfile ("config.yml"):
    config_path = "config.yml"
elif os.path.isfile ("config.json"):
    config_path = "config.json"
else:
    exit("no config file found, need config.yml or config.json")

configfile: config_path



##################################
#     DEFINE GLOBAL VARIABLES    #
##################################

# path to reference
GENOMES  = config["genome"].keys()

# raw data to analyse
DATA = config["data"]

# Define global var to know input files
INPUT_exp = list(config["data"]["input"].keys())[0]
INPUT_sample = config["data"]["input"][INPUT_exp]["fastq"]

####### variable for expand
# this variable will be used to generate target file automaticaly
# with info from config file

## peak calling
# Define kind of peak to call by experiment (broad and/or narrow)
EXP_broad = []
EXP_narrow = []

for exp in DATA["IP"].keys():
    if "narrow" in DATA["IP"][exp]["peak_calling"]["type"] :
        EXP_narrow.append(exp)
    if "broad" in DATA["IP"][exp]["peak_calling"]["type"] :
        EXP_broad.append(exp)

## peak calling for experiment with known fragment mean length
# broad and/or narrow
EXP_broad_frag = []
EXP_narrow_frag = []


for exp in DATA["IP"].keys():
    if "narrow" in DATA["IP"][exp]["peak_calling"]["type"] :
        if "params" in DATA["IP"][exp]["peak_calling"] : 
            if "fragLength" in DATA["IP"][exp]["peak_calling"]["params"] :
                EXP_narrow_frag.append(exp)
    if "broad" in DATA["IP"][exp]["peak_calling"]["type"] :
        if "params" in DATA["IP"][exp]["peak_calling"] : 
            if "fragLength" in DATA["IP"][exp]["peak_calling"]["params"] :
                EXP_broad_frag.append(exp)


# intersect peaks
# i.e.: genes.start_-1000_-500
BED_intersect_peaks   = {}

for genome in config["intersect_peaks"].keys():
    BED_intersect_peaks[genome] = []
    if genome in config["intersect_peaks"]: 
        for bedfile in config["intersect_peaks"][genome].keys():
            if bedfile in config["intersect_peaks"][genome]: 
                for ancor in config["intersect_peaks"][genome][bedfile].keys():
                    if config["intersect_peaks"][genome][bedfile][ancor] is not None :
                        for coord in config["intersect_peaks"][genome][bedfile][ancor]:
                            BED_intersect_peaks[genome].append(bedfile + "." + str(ancor) + "_" + str(coord[0]) + "_" + str(coord[1]))

######################################
# DEFINE CONFIG DEFAULT VARIABLES    #
######################################
#~ config["default_config"] if "default_config" not in config

DEFAULT_PARAMS = {
    "plot_type": {
        "upstream" : "-2000",
        "unscaled5prime" : "500",
        "regionBodyLength" : "2000",
        "unscaled3prime" : "500",
        "downstream" : "2000",
        "bin_size" : "10"
    },
    "peak_calling" :{
        "fragLength" : "200",
        "keep_dup" : "1",                # "1" or "auto" (binomial distribution) or "all"
        "call_summits" : "",             # "" or "--call-summits"
        "use_static_background" : "",    # or "--nolambda"
    },
    "bigWig":{
        "bin_size" : "10",
        "bin_smooth" : "30",
        "centerReads" : "",                         #false/--centerReads==true
        "ignoreForNormalization": "",               #chromosome to ignore
        "ignoreDuplicates" : "--ignoreDuplicates",   #true
        "minMappingQuality" : "20"
        }
}

for i in GENOMES:
    if i not in config["default_config"].keys() : 
        config["default_config"][i] = {}
    for j in DEFAULT_PARAMS.keys():
        if j not in config["default_config"][i].keys():
            config["default_config"][i][j] = {}
        for k in DEFAULT_PARAMS[j].keys():
            if k not in config["default_config"][i][j].keys():
                config["default_config"][i][j][k] = None
            if config["default_config"][i][j][k] is None:
                config["default_config"][i][j][k] = DEFAULT_PARAMS[j][k]

#~ print(config["default_config"])
#####################################
#       FUNCTIONS DEFINITION        #
#####################################
def is_empty_file(path):
    if os.path.isfile(str(path)):
        return os.stat(str(path)).st_size==0
    else:
        return 1


def simplify_list(A):
    return list(set().union(*A))

########################################################################
#### TODO ####
#
# correctGCbias
# add bigWig to Jbrowse
#
### to filter mapping Quality mapQ 20 and PCR/optical duplicate
#
# samtools view Sample_27_NoIndex.sorted.bam -b -o Sample_27_NoIndex.sorted.mapQ20.rmdup.bam -q 20 -F 1024 
#
########################################################################


rule all:
    input:
        ############################### Quality Control ###############################
        ###############################################################################
        # fastQC
        [expand("{exp}/raw_data/fastQC/{samples}.filtered.fq_fastqc.html",
                exp=EXP, samples=DATA[cond][EXP]["fastq"])
                for cond in ["input", "IP"] for EXP in DATA[cond].keys()],

        #deeptools QC
        expand("QC/{genome}/bamCorrelate.bin.spearman.svg",
                genome = GENOMES),
        expand("QC/{genome}/bamCorrelate.bin.pearson.svg",
                genome = GENOMES),
        expand("QC/{genome}/bamFingerprint.svg",
                genome = GENOMES),

        ################################### bigWig ###################################
        ##############################################################################
        #bigwig
        [expand("{exp}/bigwig/{genome}/{samples}.bw",
                exp=EXP, genome=GENOMES, samples=DATA[cond][EXP]["fastq"])
                for cond in ["input", "IP"] for EXP in DATA[cond].keys()],
        #~ #bigwig ratio over input
        [expand("{exp}/bigwig/{genome}/{samples}.input_normalised.bw",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in DATA["IP"].keys()],

        ################################ peak calling ################################
        ##############################################################################
        #narrow peak
        [expand("{exp}/peak_calling/{genome}/narrow/{samples}_{suffix}",
                exp=EXP, 
                genome=GENOMES, 
                samples=DATA["IP"][EXP]["fastq"], 
                suffix = ["peaks.narrowPeak", "control_lambda.bw", "treat_pileup.bw"])
                    for EXP in EXP_narrow],
        #broad peak
        [expand("{exp}/peak_calling/{genome}/broad/{samples}_{suffix}",
                exp=EXP, 
                genome=GENOMES, 
                samples=DATA["IP"][EXP]["fastq"], 
                suffix = ["peaks.broadPeak", "control_lambda.bw", "treat_pileup.bw"])
                    for EXP in EXP_broad],
                
        #narrow peak _with fraglen
        [expand("{exp}/peak_calling/{genome}/narrow/frag_{fragLen}/{samples}_{suffix}",
                exp=EXP, 
                genome=GENOMES, 
                samples=DATA["IP"][EXP]["fastq"], 
                fragLen = DATA["IP"][EXP]["peak_calling"]["params"]["fragLength"], 
                suffix = ["peaks.narrowPeak", "control_lambda.bw", "treat_pileup.bw"])
                    for EXP in EXP_narrow_frag],

        #broad peak _with fraglen
        [expand("{exp}/peak_calling/{genome}/broad/frag_{fragLen}/{samples}_{suffix}",
                exp=EXP, 
                genome=GENOMES, 
                samples=DATA["IP"][EXP]["fastq"], 
                fragLen = DATA["IP"][EXP]["peak_calling"]["params"]["fragLength"], 
                suffix = ["peaks.broadPeak", "control_lambda.bw", "treat_pileup.bw"])
                    for EXP in EXP_broad_frag],
                    

        ############################ heatmapper and profiler ############################
        #################################################################################
        #Heatmap and profile IP and Input(RPKM)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA[cond][EXP]["fastq"], 
                ancor = ANCOR, bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["heatmap", "profile"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for cond in ["input", "IP"] 
                                for EXP in DATA[cond].keys()],
        #Heatmap and profile IP normalized by Input (ratio over input)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.input_normalised.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA["IP"][EXP]["fastq"], 
                ancor = ANCOR, 
                bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["heatmap", "profile"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for EXP in DATA["IP"].keys()],


        ############################ peaks annotations ############################
        ###########################################################################
        #direct target
        [expand("{exp}/mat_recap/{genome}/{exp}.{samples}.{bed_intersect_peaks}.direct_target.mat_recap",
                exp=EXP, genome=GENOME, samples=DATA["IP"][EXP]["fastq"], bed_intersect_peaks = BED_intersect_peaks[GENOME])
                for GENOME in BED_intersect_peaks.keys() for EXP in DATA["IP"].keys()],

        ## TODO: re-implement GC bias
        #~ ################################ GC bias ################################
        #~ #########################################################################

        #~ #ComputeGCBias IP
        #~ [expand("{exp}/mapping/{genome}/GCBias.{samples}.png",
                #~ exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                #~ for EXP in DATA["IP"].keys()],
                
        #~ [expand("{exp}/mapping/{genome}/GCBias.{samples}.matrix.txt",
                #~ exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                #~ for EXP in DATA["IP"].keys()],

        #~ #ComputeGCBias Input
        #~ [expand("{exp}/mapping/{genome}/GCBias.Input.{samples}.png",
                #~ exp=EXP, genome=GENOMES, samples=DATA["input"][EXP]["fastq"])
                #~ for EXP in DATA["input"].keys()],
        #~ [expand("{exp}/mapping/{genome}/GCBias.Input.{samples}.matrix.txt",
                #~ exp=EXP, genome=GENOMES, samples=DATA["input"][EXP]["fastq"])
                #~ for EXP in DATA["input"].keys()],





# target rules, part of rule all
rule QC:
    input:
        ############################### Quality Control ###############################
        ###############################################################################
        # fastQC
        [expand("{exp}/raw_data/fastQC/{samples}.filtered.fq_fastqc.html",
                exp=EXP, samples=DATA[cond][EXP]["fastq"])
                for cond in ["input", "IP"] for EXP in DATA[cond].keys()],

        #deeptools QC
        expand("QC/{genome}/bamCorrelate.bin.spearman.svg",
                genome = GENOMES),
        expand("QC/{genome}/bamCorrelate.bin.pearson.svg",
                genome = GENOMES),
        expand("QC/{genome}/bamFingerprint.svg",
                genome = GENOMES)
rule heatmap:
    input:
        ############################ heatmapper ############################
        ####################################################################
        #Heatmap and profile IP and Input(RPKM)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA[cond][EXP]["fastq"], 
                ancor = ANCOR, bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["heatmap"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for cond in ["input", "IP"] 
                                for EXP in DATA[cond].keys()],
        #Heatmap and profile IP normalized by Input (ratio over input)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.input_normalised.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA["IP"][EXP]["fastq"], 
                ancor = ANCOR, 
                bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["heatmap"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for EXP in DATA["IP"].keys()]
rule profiles:
    input:
        ############################ profiler ############################
        ##################################################################
        #IP and Input(RPKM)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA[cond][EXP]["fastq"], 
                ancor = ANCOR, bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["profile"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for cond in ["input", "IP"] 
                                for EXP in DATA[cond].keys()],
        #IP normalized by Input (ratio over input)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.input_normalised.{bed}.{ancor}.{plot_type}.png",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA["IP"][EXP]["fastq"], 
                ancor = ANCOR, 
                bed = config["plot_type"][GENOME][ANCOR], 
                plot_type = ["profile"])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for EXP in DATA["IP"].keys()],
rule mapping:
    input:
        #################################### mapping ###################################
        ################################################################################
        [expand("{exp}/mapping/{genome}/{samples}.sorted.{bam}",
                exp=EXP, genome=GENOMES, samples=DATA[cond][EXP]["fastq"], bam = ["bam", "bam.bai"])
                for cond in ["input", "IP"] for EXP in DATA[cond].keys()],
rule bigwig:
    input:
        ################################### bigWig ###################################
        ##############################################################################
        #bigwig
        [expand("{exp}/bigwig/{genome}/{samples}.bw",
                exp=EXP, genome=GENOMES, samples=DATA[cond][EXP]["fastq"])
                for cond in ["input", "IP"] for EXP in DATA[cond].keys()],
        #~ #bigwig ratio over input
        [expand("{exp}/bigwig/{genome}/{samples}.input_normalised.bw",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in DATA["IP"].keys()],
rule peak_calling:
    input:
        ################################ peak calling ################################
        ##############################################################################
        #narrow peak
        [expand("{exp}/peak_calling/{genome}/narrow/{samples}_peaks.narrowPeak",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in EXP_narrow],

        #broad peak
        [expand("{exp}/peak_calling/{genome}/broad/{samples}_peaks.broadPeak",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in EXP_broad],
                
        #narrow peak _with fraglen
        [expand("{exp}/peak_calling/{genome}/narrow/frag_{fragLen}/{samples}_peaks.narrowPeak",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"], fragLen = DATA["IP"][EXP]["peak_calling"]["params"]["fragLength"])
                for EXP in EXP_narrow_frag],

        #broad peak _with fraglen
        [expand("{exp}/peak_calling/{genome}/broad/frag_{fragLen}/{samples}_peaks.broadPeak",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"], fragLen = DATA["IP"][EXP]["peak_calling"]["params"]["fragLength"])
                for EXP in EXP_broad_frag]



# target rules, not in rule all
rule RDS:
    input:
        ############################ heatmapper and profiler ############################
        #################################################################################
        #Heatmap and profile IP and Input(RPKM)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.{bed}.{ancor}.RDS",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA[cond][EXP]["fastq"], 
                ancor = ANCOR,
                bed = config["plot_type"][GENOME][ANCOR])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for cond in ["input", "IP"] 
                                for EXP in DATA[cond].keys()],
        #Heatmap and profile IP normalized by Input (ratio over input)
        [expand("{exp}/profile/{genome}/{bed}/{ancor}/{exp}.{samples}.input_normalised.{bed}.{ancor}.RDS",
                exp=EXP, 
                genome=GENOME, 
                samples=DATA["IP"][EXP]["fastq"], 
                ancor = ANCOR, 
                bed = config["plot_type"][GENOME][ANCOR])
                    for GENOME in config["plot_type"].keys() 
                        for ANCOR in config["plot_type"][GENOME].keys() 
                            for EXP in DATA["IP"].keys()]
rule xcor:
    input:
        ################################ SPP xcor ################################
        ##########################################################################
        [expand("{exp}/peak_calling/{genome}/SPP_xcor/{samples}.spp.xcor.pdf",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in DATA["IP"]],

        ################################ MaSC xcor ###############################
        ##########################################################################
        [expand("{exp}/peak_calling/{genome}/MaSC_xcor/{samples}.xcor_MaSC.png",
                exp=EXP, genome=GENOMES, samples=DATA["IP"][EXP]["fastq"])
                for EXP in DATA["IP"]],
rule direct_target:
    input:
        # direct target
        [expand("{exp}/mat_recap/{genome}/{exp}.{samples}.{bed_intersect_peaks}.direct_target.mat_recap",
            exp=EXP, genome=GENOME, samples=DATA["IP"][EXP]["fastq"], bed_intersect_peaks = BED_intersect_peaks[GENOME])
            for GENOME in BED_intersect_peaks.keys() for EXP in DATA["IP"].keys()],
        # merge_bed
        [expand("{exp}/peak_calling/{genome}/{samples}_peaks.merged.bed",
            exp=EXP, genome=GENOME, samples=DATA["IP"][EXP]["fastq"])
            for GENOME in BED_intersect_peaks.keys() for EXP in DATA["IP"].keys()]

###############################
#       QUALITY CONTROL       #
###############################

###### Sequencing Quality ######
################################

rule fast_QC:
    input:
        "{exp}/raw_data/{samples}.filtered.fq.gz"
    output:
        "{exp}/raw_data/fastQC/{samples}.filtered.fq_fastqc.html",
        "{exp}/raw_data/fastQC/{samples}.filtered.fq_fastqc.zip"
    threads: 28
    log:
        "{exp}/raw_data/fastQC/{samples}.fastqc.log"
    run:
        shell("mkdir -p {wildcards.exp}/raw_data/fastQC/")
        shell("fastqc -o {wildcards.exp}/raw_data/fastQC/ -t {threads} {input} 2>{log}")


## strand cross correlation ##
##############################
#
# Determine mean fragments length
# and quality of IP
#

# xcor with SPP
rule spp_xcor:
    input:  bam = "{exp}/mapping/{genome}/{samples}.sorted.bam",
            bai = "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output: "{exp}/peak_calling/{genome}/SPP_xcor/{samples}.spp.xcor.pdf"
    threads: 28
    params:
        maxShift = 300,
        bin = 1
    log : "{exp}/peak_calling/{genome}/SPP_xcor/{samples}.spp.xcor.log"
    shell : "Rscript " + config["SPP"]["xcor"] + " {input.bam} {output} {threads} {params.maxShift} {params.bin} >{log}"

# xcor with MaSC (taking acount of mappability)
rule MaSC_xcor:
    input:  bam = "{exp}/mapping/{genome}/{samples}.sorted.bam",
            bai = "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output: "{exp}/peak_calling/{genome}/MaSC_xcor/{samples}.xcor_MaSC.png"
    threads: 1
    params:
        mappability_path = lambda wildcards : config["genome"][wildcards.genome]["mappability"],
        chrom_length_file = lambda wildcards : config["genome"][wildcards.genome]["chr_len"],
        smooth_win_size = 15,
        min_shift = 0,
        max_shift = 400,
    log : "{exp}/peak_calling/{genome}/MaSC_xcor/{samples}.xcor_MaSC.log"
    shell : "bedtools bamtobed -i {input.bam} | "
            "MaSC.pl --verbose "
            "--mappability_path {params.mappability_path} "
            "--chrom_length_file {params.chrom_length_file} "
            "--prefix {wildcards.exp}/peak_calling/{wildcards.genome}/MaSC_xcor/{wildcards.samples}.xcor "
            "--smooth_win_size {params.smooth_win_size} "
            "--min_shift {params.min_shift} "
            "--max_shift {params.max_shift} "
            "--input_bed - "
            ">{log}"

# Know if IP is broad or narrow and input well sequenced
rule bamFingerprint:
    input:
        bamFiles = lambda wildcards: simplify_list([expand(
                    "{exp}/mapping/{genome}/{samples}.sorted.bam",
                    exp=EXP, genome=wildcards.genome, samples=DATA[cond][EXP]["fastq"])
                    for cond in ["input", "IP"] for EXP in DATA[cond].keys()]),
        baiFiles = lambda wildcards: simplify_list([expand(
                    "{exp}/mapping/{genome}/{samples}.sorted.bam.bai",
                    exp=EXP, genome=wildcards.genome, samples=DATA[cond][EXP]["fastq"])
                    for cond in ["input", "IP"] for EXP in DATA[cond].keys()])
    output:
        "QC/{genome}/bamFingerprint.svg"
    priority: 
        5
    params:
        ignoreDuplicates = "",
        minMappingQuality = "20",
        binSize = "500",
        fragmentLength = lambda wildcards: config["default_config"][wildcards.genome]["peak_calling"]["fragLength"],
        numberOfSamples = "500000", # Number of bins that sampled from the genome,
        plotTitle = "'Fingerprint of ChIP-seq data'"
    threads: 28
    run:
        tab = []
        # get the name of first directory
        for i in input.bamFiles:
            m = re.search('(^[^/]+)/', i)
            tab.append(m.group(1))
        
        tab = " ".join(tab)
        
        shell(
            "plotFingerprint --bamfiles {input.bamFiles} "
            "--plotFile {output} "
            "--labels "+tab+" "
            "--numberOfProcessors {threads} "
            "{params.ignoreDuplicates} "
            "--minMappingQuality {params.minMappingQuality} "
            "--binSize {params.binSize} "
#            "--fragmentLength {params.fragmentLength} " # deeptools  version 1.5
            "--extendReads {params.fragmentLength} " # deeptools version 2.0
            "--numberOfSamples {params.numberOfSamples} "
            "--plotTitle {params.plotTitle} "
        )

####### bam correlate #######
#############################
#
# correlation between replicate and vs input/other IP
# 

rule bamCorrelate:
    input:
        bamFiles = lambda wildcards: simplify_list([expand(
                    "{exp}/mapping/{genome}/{samples}.sorted.bam",
                    exp=EXP, genome=wildcards.genome, samples=DATA[cond][EXP]["fastq"])
                    for cond in ["input", "IP"] for EXP in DATA[cond].keys()]),
        baiFiles = lambda wildcards: simplify_list([expand(
                    "{exp}/mapping/{genome}/{samples}.sorted.bam.bai",
                    exp=EXP, genome=wildcards.genome, samples=DATA[cond][EXP]["fastq"])
                    for cond in ["input", "IP"] for EXP in DATA[cond].keys()]),
        exp = lambda wildcards: simplify_list([expand(
                    "{exp}",
                    exp=EXP)
                    for cond in ["input", "IP"] for EXP in DATA[cond].keys()])
    output:
        "QC/{genome}/bamCorrelate.bin.npz"
    threads: 28
    version: subprocess.getoutput("multiBamSummary --version")
    log:
        "QC/{genome}/bamCorrelate.log"
    run:
        tab = []
        for i in input.bamFiles:
            m = re.search('(^[^/]+)/', i)
            tab.append(m.group(1))
        
        tab = " ".join(tab)        
        shell("multiBamSummary bins --bamfiles {input.bamFiles} --labels "+tab+" --numberOfProcessors {threads} -out {output} 2>{log}")

rule plotCorrelation:
    input:
        "QC/{genome}/bamCorrelate.bin.npz"
    output:
        "QC/{genome}/bamCorrelate.bin.{corMethod}.svg"
    version: subprocess.getoutput("plotCorrelation --version")
    log:
        "QC/{genome}/bamCorrelate.bin.{corMethod}.log"
    threads:
        1
    params:
        whatToPlot = "heatmap",#{heatmap,scatterplot}
        skipZeros = "", #--skipZeros
        plotTitle = lambda wildcards: "'"+wildcards.corMethod+" correlation'",
        removeOutliers = "",#--removeOutliers
        plotNumbers = "--plotNumbers" # if whatToPlot = "heatmap"
    shell:
        "plotCorrelation -in {input} "
        "-o {output} "
        "--whatToPlot {params.whatToPlot} "
        "{params.skipZeros} "
        "{params.removeOutliers} "
        "{params.plotNumbers} "
        "--corMethod {wildcards.corMethod} "
        "--plotTitle {params.plotTitle} "
        "&> {log}"

#########################################################
#########################################################
#            # # # # # # # # # # # # # # # #            #
#            #       Data processing       #            #
#            # # # # # # # # # # # # # # # #            #
#########################################################
#########################################################

# /!\ TO TEST: if file name is not in good format generate sym link 
rule fastq_gz_to_fq_gz_extension:
    input:"{base}.fastq.gz"
    output:"{base}.fq.gz"
    threads: 1
    shell:"ln -s {input} {output}"

rule gzip:
    input:"{base}"
    output:"{base}.gz"
    threads: 1
    shell:"gzip {input}"


###############################
#  Filter good quality reads  #
###############################

rule illumina_filtering:
    input: "{samples}.fq.gz"
    output: "{samples}.filtered.fq.gz"
    log: "{samples}.fastq_filtering.log"
    threads: 3
    shell:  "gunzip -c {input} | "
            "fastq_illumina_filter -vvN 2>{log} | "
            "gzip > {output}"


###############################
#         bwa mapping         #
#     sort/index BAM files    #
###############################
rule bwa_index:
    input:
        "{genome}"
    output:
        "{genome}.bwt"
    message:
        "bwa index of {wildcards.genome}"
    threads: 1
    shell:
        "bwa index {input}"

rule bwa_aln:
    input:
        bwt = lambda wildcards: config["genome"][wildcards.genome]["index"] + '.bwt',
        fasta = lambda wildcards: config["genome"][wildcards.genome]["index"],
        fq = "{exp}/raw_data/{samples}.filtered.fq.gz"
    output:
        temp("{exp}/mapping/{genome, [^/]+}/{samples}.sai")
    version: subprocess.getoutput("bwa 2>&1|grep Version:|sed -r 's/Version: +//'")
    threads: 20
    log:
        "{exp}/mapping/{genome}/{samples}.bwa_aln.log"
    shell:
        "bwa aln -t {threads} {input.fasta} {input.fq} >{output} 2>{log}"

rule bwa_sampe_to_bam:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["index"],
        sai1 = "{exp}/mapping/{genome}/{samples}_1.sai",
        fq1 = "{exp}/raw_data/{samples}_1.filtered.fq.gz",
        sai2 = "{exp}/mapping/{genome}/{samples}_2.sai",
        fq2 = "{exp}/raw_data/{samples}_2.filtered.fq.gz"
    output:
        temp("{exp}/mapping/{genome, [^/]+}/{samples}.bam")
    log:
        "{exp}/mapping/{genome}/{samples}.samse.log"
    threads: 2
    shell:"bwa sampea {input.fasta} {input.sai1} {input.sai2} {input.fq1} {input.fq2} 2>{log}| samtools view -Sbh - > {output}"

rule bwa_samse_to_bam:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["index"],
        sai = "{exp}/mapping/{genome}/{samples}.sai",
        fq = "{exp}/raw_data/{samples}.filtered.fq.gz"
    output:
        temp("{exp}/mapping/{genome}/{samples}.bam")
    log:
        "{exp}/mapping/{genome}/{samples}.samse.log"
    threads: 2
    shell:"bwa samse {input.fasta} {input.sai} {input.fq} 2>{log}| samtools view -Sbh - > {output}"

rule sort_bam:
    input:
        "{base}.bam"
    output:
        "{base}.sorted.bam"
    threads: 10
    log:
        "{base}.sort.log"
    shell:
        "samtools sort -l 9 -m 2G -@ {threads} {input} -o {output} 2>{log}"

rule index_bam:
    input:
        "{base}.sorted.bam"
    output:
        "{base}.sorted.bam.bai"
    priority: 50
    threads: 1
    shell:
        "samtools index {input}"

###############################
#        BigWig MAKING        #
###############################
rule input_normalisation_bw:
    input:
        "{exp}/bigwig/{genome}/{samples}.bw",
        INPUT_exp + "/bigwig/{genome}/" + INPUT_sample + ".bw"
    output:
        "{exp}/bigwig/{genome}/{samples}.input_normalised.bw"
    threads: 28
    version: subprocess.getoutput("bigwigCompare --version")
    params:
        pseudocount = "1",
        ratio = "ratio", #{log2,ratio,subtract,add,reciprocal_ratio}
        binSize = "10"
    run:
        if input[1] != input[0]:
            shell(
                "bigwigCompare "
                "--bigwig1 {input[0]} "
                "--bigwig2 {input[1]} "
                "--outFileName {output} "
                "--numberOfProcessors {threads} "
                "--pseudocount {params.pseudocount} "
                "--ratio {params.ratio} "
                "--binSize {params.binSize} "
            )


rule bam_to_bigWig:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output:
        "{exp}/bigwig/{genome}/{samples}.bw"
    log:
        "{exp}/bigwig/{genome}/{samples}.bw.log"
    threads: 28
    version: subprocess.getoutput("bamCoverage --version")
    params:
        binSize = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["bin_size"],
        smoothLength = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["bin_smooth"],
        centerReads = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["centerReads"],
        ignoreForNormalization = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["ignoreForNormalization"],
        ignoreDuplicates = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["ignoreDuplicates"],
        minMappingQuality = lambda wildcards: config["default_config"][wildcards.genome]["bigWig"]["minMappingQuality"]
    run:
        extendReads = config["default_config"][wildcards.genome]["peak_calling"]["fragLength"] # default value
        
        if (wildcards.exp in DATA["IP"] and
            "peak_calling" in DATA["IP"][wildcards.exp] and
            "params" in DATA["IP"][wildcards.exp]["peak_calling"] and 
            "fragLength" in DATA["IP"][wildcards.exp]["peak_calling"]["params"]): 
            extendReads = DATA["IP"][wildcards.exp]["peak_calling"]["params"]["fragLength"]
        elif (wildcards.exp in DATA["input"] and 
            "peak_calling" in DATA["input"][wildcards.exp] and
            "params" in DATA["input"][wildcards.exp]["peak_calling"] and 
            "fragLength" in DATA["input"][wildcards.exp]["peak_calling"]["params"]):
            extendReads = DATA["input"][wildcards.exp]["peak_calling"]["params"]["fragLength"]
        shell(
            "bamCoverage "
            "-b {input[0]} "
            "-o {output} "
            "--normalizeUsingRPKM "
            "--numberOfProcessors {threads} "
            "--binSize {params.binSize} "
            "{params.centerReads} "
            "{params.ignoreForNormalization} "
            "--extendReads " + str(extendReads) + " "     # for version 2.0
#            "--fragmentLength " + extendReads + " "   # for version 1.0
            "--smoothLength {params.smoothLength} "
            "{params.ignoreDuplicates} "
            "--minMappingQuality {params.minMappingQuality} "
            "--verbose &>{log}"
        )


###############################
#        PEAKS CALLING        #
###############################

rule macs2_peak_calling_frag_length:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        INPUT_exp + "/mapping/{genome}/" + INPUT_sample + ".sorted.bam"
    output:
        "{exp}/peak_calling/{genome}/{peak_type}/frag_{fragLen}/{samples}_peaks.{peak_type}Peak",
        temp("{exp}/peak_calling/{genome}/{peak_type}/frag_{fragLen}/{samples}_control_lambda.bdg"), 
        temp("{exp}/peak_calling/{genome}/{peak_type}/frag_{fragLen}/{samples}_treat_pileup.bdg")
    log:
        "{exp}/peak_calling/{genome}/{peak_type}/frag_{fragLen}/{samples}_macs2.{peak_type}.log"
    # use default parameter or parameter for specific IP (depending on config.yaml)
    params:
        peak_type = lambda wildcards : 
            "--broad" if wildcards.peak_type == "broad" 
            else "",
        genome_size = lambda wildcards :
            config["genome"][wildcards.genome]["effective_size"],
        use_static_background =  lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["use_static_background"] # "[]" or "--nolambda"
                if "peak_calling" in config["data"]["IP"][wildcards.exp] and
                    "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "use_static_background" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["use_static_background"],
        call_summits = lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["call_summits"] # "[]" or "--call-summits"
                if "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "call_summits" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["call_summits"],
        keep_dup = lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["keep_dup"] # "[1]" or "auto" (binomial distribution) or "all"
                if "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "keep_dup" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["keep_dup"],
        nomodel = "--nomodel",
        extsize = lambda wildcards: wildcards.fragLen,
        
    threads: 1
    run:
        if input[1] != input[0]:
            shell(
                "macs2 callpeak "
                "-t {input[0]} "
                "-c {input[1]} "
                "-f BAM "
		"--tempdir /home/jean-philippe.villemin/tmp_deeptools/ "
                "--name {wildcards.samples} "
                "--gsize {params.genome_size} "
                "{params.nomodel} " 
                "--extsize {params.extsize} "
                "{params.peak_type} "
                "{params.use_static_background} "
                "{params.call_summits} "
                "--keep-dup {params.keep_dup} "
                "--bdg "
                "--outdir {wildcards.exp}/peak_calling/{wildcards.genome}/{wildcards.peak_type}/frag_{wildcards.fragLen} "
                "2>{log}"
            )
            ## verify if macs have generated files (at least 1 peak called)
            ## if not create empty file to avoid snakemake stop
            if is_empty_file(wildcards.exp + "/peak_calling/" + wildcards.genome + "/" + wildcards.peak_type + "/" + wildcards.samples + "_peaks." + wildcards.peak_type + "Peak"):
                shell(
                    "touch {wildcards.exp}/peak_calling/{wildcards.genome}/{wildcards.peak_type}/frag_{wildcards.fragLen}/{wildcards.samples}_peaks.{wildcards.peak_type}Peak"
                )

rule macs2_peak_calling_no_frag_length:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        INPUT_exp + "/mapping/{genome}/" + INPUT_sample + ".sorted.bam"
    output:
        "{exp}/peak_calling/{genome}/{peak_type}/{samples}_peaks.{peak_type}Peak",
        temp("{exp}/peak_calling/{genome}/{peak_type}/{samples}_control_lambda.bdg"), 
        temp("{exp}/peak_calling/{genome}/{peak_type}/{samples}_treat_pileup.bdg")
    log:
        "{exp}/peak_calling/{genome}/{peak_type}/{samples}_macs2.{peak_type}.log"
    # use default parameter or parameter for specific IP (depending on config.yaml)
    params:
        peak_type = lambda wildcards : 
            "--broad" if wildcards.peak_type == "broad" 
            else "",
        genome_size = lambda wildcards :
            config["genome"][wildcards.genome]["effective_size"],
        use_static_background =  lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["use_static_background"] # "[]" or "--nolambda"
                if "peak_calling" in config["data"]["IP"][wildcards.exp] and
                    "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "use_static_background" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["use_static_background"],
        call_summits = lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["call_summits"] # "[]" or "--call-summits"
                if "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "call_summits" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["call_summits"],
        keep_dup = lambda wildcards :
            config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]["keep_dup"] # "[1]" or "auto" (binomial distribution) or "all"
                if "params" in config["data"]["IP"][wildcards.exp]["peak_calling"] and
                    "keep_dup" in config["data"]["IP"][wildcards.exp]["peak_calling"]["params"]
                else
            config["default_config"][wildcards.genome]["peak_calling"]["keep_dup"],
        nomodel = ""        
    threads: 1
    run:
        if input[1] != input[0]:
            shell(
                "macs2 callpeak "
                "-t {input[0]} "
                "-c {input[1]} "
                "-f BAM "
		"--tempdir /home/jean-philippe.villemin/tmp_deeptools/ "
                "--name {wildcards.samples} "
                "--gsize {params.genome_size} "
                "{params.nomodel} " 
                "{params.peak_type} "
                "{params.use_static_background} "
                "{params.call_summits} "
                "--keep-dup {params.keep_dup} "
                "--bdg "
                "--outdir {wildcards.exp}/peak_calling/{wildcards.genome}/{wildcards.peak_type} "
                "2>{log}"
            )
            ## verify if macs have generated files (at least 1 peak called)
            ## if not create empty file to avoid snakemake stop
            if is_empty_file(wildcards.exp + "/peak_calling/" + wildcards.genome + "/" + wildcards.peak_type + "/" + wildcards.samples + "_peaks." + wildcards.peak_type + "Peak"):
                shell(
                    "touch {wildcards.exp}/peak_calling/{wildcards.genome}/{wildcards.peak_type}/{wildcards.samples}_peaks.{wildcards.peak_type}Peak"
                )

rule macs2_bedGraph_to_bigwig:
    input:"{exp}/peak_calling/{genome}/{peak_type}/{suffix}.bdg"
    output:"{exp}/peak_calling/{genome}/{peak_type, (narrow)|(broad)}/{suffix}.bw"
    params:
        chrom_size = lambda wildcards: config["genome"][wildcards.genome]["chr_len"]
    threads:1
    log:
    shell:"bedGraphToBigWig {input} {params.chrom_size} {output}"

## merge all peaks called ##
ruleorder: 
    merge_broad_narow_peaks > 
    make_symbolic_link_narrow_peaks > 
    make_symbolic_link_broad_peaks

rule merge_broad_narow_peaks:
    input:
        "{exp}/peak_calling/{genome}/broad/{samples}_peaks.broadPeak",
        "{exp}/peak_calling/{genome}/narrow/{samples}_peaks.narrowPeak"
    output:
        "{exp}/peak_calling/{genome}/{samples}_peaks.merged.bed"
    threads : 2
    shell:
        "cat {wildcards.exp}/peak_calling/{wildcards.genome}/*/{wildcards.samples}_peaks.[nb]*Peak "
        "{wildcards.exp}/peak_calling/{wildcards.genome}/*/*/{wildcards.samples}_peaks.[nb]*Peak | "
        "cut -f1-3 | "
        "bedtools sort -i - | "
        "bedtools merge -i - >{output}"

rule make_symbolic_link_broad_peaks:
    input:
        "{exp}/peak_calling/{genome}/broad/{samples}_peaks.broadPeak"
    output:
        "{exp}/peak_calling/{genome}/{samples}_peaks.merged.bed"
    threads : 1
    shell:
        "ls -s {input} {output}"

rule make_symbolic_link_narrow_peaks:
    input:
        "{exp}/peak_calling/{genome}/narrow/{samples}_peaks.narrowPeak"
    output:
        "{exp}/peak_calling/{genome}/{samples}_peaks.merged.bed"
    threads : 1
    shell:
        "ls -s {input} {output}"


#### Make mat_recap files   ####
################################
#
# Make matrix of boolean to know if features in bed have a peak
# in or around it's coordinates

rule direct_target_peaks:
    input:
        "{exp}/peak_calling/{genome}/{samples}_peaks.merged.bed"
    output:
        "{exp}/mat_recap/{genome}/{exp}.{samples}.{bed}.{intersection}.direct_target.mat_recap"
    threads : 1
    run:
        bed_key = wildcards.bed
        # récupère les paramètres pour l'intersection
        tab = wildcards.intersection.split('_')
        if len(tab) == 3:
            ancor,up,down = tab
        else:
            print("error: verify your config file (intersect_peaks parameters)")
            exit
        
        #crée le dossier contenant les fichiers mat_racap
        shell("mkdir -p {wildcards.exp}/mat_recap/{wildcards.genome}")
        
        
        # choisi la commande en fonction de l'ancre (start/TSS, end/TES, body, center))
        # si start (2e col) devient negatif, alors start = 0 (pour éviter des erreur avec bedtools)
        if ancor in ["start", "TSS"]:
            awk_cmd = 'awk \'OFS="\t" {{if($6 == "-"){{TSS = $3; $2 = TSS - 1 - '+down+'; $3 = TSS - '+up+'}}else{{TSS = $2; $2 = TSS + '+up+' ;$3 = TSS + 1 + '+down+'}}if($2<0){{$2=0}}print}}\' '
        elif ancor in ["end", "TES"]:
            awk_cmd = 'awk \'OFS="\t" {{if($6 == "-"){{TES = $2; $2 = TES - '+down+'; $3 = TES + 1 - '+up+'}}else{{TES = $3; $2 = TES - 1 + '+up+' ;$3 = TES + '+down+'}}if($2<0){{$2=0}}print}}\' '
        elif ancor == "center":
            awk_cmd = 'awk \'OFS="\t" {{if($6 == "-"){{center = int(($3+$2)/2); $2 = center - '+down+'; $3 = center + 1 - '+up+'}}else{{center = int(($3+$2)/2); $2 = center + '+up+' ;$3 = center + 1 + '+down+'}}if($2<0){{$2=0}}print}}\' '
        elif ancor == "body":
            awk_cmd = 'awk \'OFS="\t" {{if($6 == "-"){{$2 = $2 - '+down+'; $3 = $3 - '+up+'}}else{{$2 = $2 + '+up+' ;$3 = $3 + '+down+'}}if($2<0){{$2=0}}print}}\' '

        # suite de la commande
        bedtool_cmd = "bedtools intersect -a - -b {input} -wa | cut -f4 > " + wildcards.exp + "/mat_recap/" + wildcards.genome + "/" + wildcards.bed + "."+ wildcards.intersection +".positive_features.tmp"

        # création du fichier temporaire contenant les features positives
        shell(awk_cmd + config["bed_files"][bed_key] + " | " + bedtool_cmd)

        # on crée un objet set contenant toutes les features positives
        positive_features = set(line.strip() for line in open(wildcards.exp+"/mat_recap/"+wildcards.genome+"/"+wildcards.bed+"."+wildcards.intersection+".positive_features.tmp", "r"))
        
        ##make pseudo mat_recap
        #prepare file
        shell(
            "cut -f4 "+config["bed_files"][bed_key]+" >{wildcards.exp}/mat_recap/{wildcards.genome}/{wildcards.bed}.{wildcards.intersection}.acc.tmp"
        )

        bed_file = open(wildcards.exp+"/mat_recap/"+wildcards.genome+"/"+wildcards.bed+"."+wildcards.intersection+".acc.tmp", "r")
        mat_recap = open(output[0], "w")


        ## make mat_recap
        mat_recap.write("ft_acc\t"+wildcards.bed+"."+wildcards.intersection+"\n")

        for feature in bed_file:
            feature = feature.strip()
            line = feature + "\t"
            # if  is a positive feature then add 1 else add 0
            line = line + str(int(feature in positive_features)) + "\n"
            mat_recap.write(line)

        # rm tempory files
        shell(
        "rm {wildcards.exp}/mat_recap/{wildcards.genome}/{wildcards.bed}.{wildcards.intersection}*tmp"
        )


#########################################################
#########################################################
#            # # # # # # # # # # # # # # # #            #
#            #     data visualisation      #            #
#            # # # # # # # # # # # # # # # #            #
#########################################################
#########################################################


###############################
#        computeMatrix        #
###############################

rule computeMatrix_ref_point:
    input:
        "{exp}/bigwig/{genome}/{samples}.bw"
    output:
        "{exp}/profile/{genome}/{bed}/{ref, (TSS)|(TES)|(center)}/{exp}.{samples}.{bed}.{ref}.matrix.gz"
    version: subprocess.getoutput("computeMatrix --version")
    params:
        upstream = lambda wildcards : 
                    config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]["upstream"]
                        if "upstream" in config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["upstream"],
        downstream = lambda wildcards : 
                    config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]["downstream"] 
                        if "downstream" in config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["downstream"],
        binSize = lambda wildcards : 
                    config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]["bin_size"] 
                        if "bin_size" in config["plot_type"][wildcards.genome][wildcards.ref][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["bin_size"],
        sortRegions = "descend",#{descend,ascend,no}
        sortUsing = "mean"
    threads: 28
    log:
        "{exp}/profile/{genome}/{bed}/{ref}/{exp}.{samples}.{bed}.{ref}.matrix.log"
    run:
        bed_key = wildcards.bed
        referencePoint = wildcards.ref
        shell(
            "computeMatrix reference-point "
            "--regionsFileName " + config["bed_files"][bed_key] + " "
            "--scoreFileName {input} "
            "--outFileName {output} "
            "--referencePoint "+ referencePoint +" "
            "--upstream {params.upstream} "
            "--downstream {params.downstream} "
            "--binSize {params.binSize} "
            "--sortRegions {params.sortRegions} "
            "--sortUsing {params.sortUsing} "
            "--numberOfProcessors {threads} "
            "2>{log} "
        )

rule computeMatrix_unscalled_body:
    input:
        "{exp}/bigwig/{genome}/{samples}.bw"
    output:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.matrix.gz"
    version: subprocess.getoutput("computeMatrix --version")
    params:
        upstream = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["upstream"]
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "upstream" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["upstream"],
        unscaled5prime = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["unscaled5prime"] 
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "unscaled5prime" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["unscaled5prime"],
        regionBodyLength = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["regionBodyLength"] 
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "regionBodyLength" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["regionBodyLength"],
        unscaled3prime = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["unscaled3prime"] 
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "unscaled3prime" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["unscaled3prime"],
        downstream = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["downstream"] 
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "downstream" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["downstream"],
        binSize = lambda wildcards : 
                    config["plot_type"][wildcards.genome]["body"][wildcards.bed]["bin_size"] 
                        if config["plot_type"][wildcards.genome]["body"][wildcards.bed] is not None and
                        "bin_size" in config["plot_type"][wildcards.genome]["body"][wildcards.bed]
                        else
                    config["default_config"][wildcards.genome]["plot_type"]["bin_size"],
        sortRegions = "descend",#{descend,ascend,no}
        sortUsing = "mean"
    threads: 28
    log:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.matrix.log"
    run:
        bed_key = wildcards.bed
        shell(
            "computeMatrix scale-regions "
            "--regionsFileName " + config["bed_files"][bed_key] + " "
            "--scoreFileName {input} "
            "--outFileName {output} "
            "--upstream {params.upstream} "
            "--unscaled5prime {params.unscaled5prime} "
            "--regionBodyLength {params.regionBodyLength} "
            "--unscaled3prime {params.unscaled3prime} "            
            "--downstream {params.downstream} "
            "--binSize {params.binSize} "
            "--sortRegions {params.sortRegions} "
            "--sortUsing {params.sortUsing} "
            "--numberOfProcessors {threads} "
            "2>{log} "
        )

###############################
#    convert matrix to RDS    #
###############################

rule deeptools_Matrix_to_RDS:
    input:
        "{base}.matrix.gz"
    output:
        "{base}.RDS"
    version: "1.0"
    params:
    threads: 1
    log:
        "{base}.deeptools_Matrix_to_RDS.log"
    shell:
        "Rscript "
        + config["hm_RDS"] + " "
        "{input} "
        "{output} "
        ">{log} "

###############################
#         heatmapper          #
###############################
rule headmapper_ref_point:
    input:
        "{exp}/profile/{genome}/{bed}/{ref, (TSS)|(TES)|(center)}/{exp}.{samples}.{bed}.{ref}.matrix.gz"
    output:
        "{exp}/profile/{genome}/{bed}/{ref, (TSS)|(TES)|(center)}/{exp}.{samples}.{bed}.{ref}.heatmap.png"
    params:
        kmeans = "", #"--kmeans 8",
        sortRegions = "descend", #{descend,ascend,no}
        sortUsing = "mean", #{mean,median,max,min,sum,region_length}
        averageTypeSummaryPlot = "mean", #{mean,median,min,max,std,sum}
        zMin = "", #"--zMin 0"
        zMax = "", #"--zMax 100"
        heatmapHeight = "25",
        heatmapWidth = "7.5",
        whatToShow = "'heatmap and colorbar'" # {plot, heatmap and colorbar,plot and heatmap,heatmap only,colorbar only,heatmap and colorbar}
    threads: 1
    run:
        refPointLabel = wildcards.ref
        bedFile = wildcards.bed
        shell(
        "plotHeatmap --matrixFile {input} --outFileName {output}"
        "{params.kmeans} "
        "--sortRegions {params.sortRegions} "
        "--sortUsing {params.sortUsing} "
        "--averageTypeSummaryPlot {params.averageTypeSummaryPlot} "
        "{params.zMin} "
        "{params.zMax} "
        "--heatmapHeight {params.heatmapHeight} "
        "--heatmapWidth {params.heatmapWidth} "
        "--whatToShow {params.whatToShow} "
        "--refPointLabel {wildcards.ref} "
        "--plotTitle '{wildcards.bed} around {wildcards.ref}: {wildcards.exp}' "
        )

rule headmapper_body:
    input:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.matrix.gz"
    output:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.heatmap.png"
    params:
        kmeans = "", #"--kmeans 8",
        sortRegions = "descend", #{descend,ascend,no}
        sortUsing = "mean", #{mean,median,max,min,sum,region_length}
        averageTypeSummaryPlot = "mean", #{mean,median,min,max,std,sum}
        zMin = "", #"--zMin 0"
        zMax = "", #"--zMax 100"
        heatmapHeight = "25",
        heatmapWidth = "7.5",
        whatToShow = "'heatmap and colorbar'" # {plot, heatmap and colorbar,plot and heatmap,heatmap only,colorbar only,heatmap and colorbar}
    threads: 1
    shell:
        "plotHeatmap --matrixFile {input} --outFileName {output} "
        "{params.kmeans} "
        "--sortRegions {params.sortRegions} "
        "--sortUsing {params.sortUsing} "
        "--averageTypeSummaryPlot {params.averageTypeSummaryPlot} "
        "{params.zMin} "
        "{params.zMax} "
        "--heatmapHeight {params.heatmapHeight} "
        "--heatmapWidth {params.heatmapWidth} "
        "--whatToShow {params.whatToShow} "
        "--plotTitle 'body of {wildcards.bed} on {wildcards.exp}' "


###############################
#          profiler           #
###############################
rule profiler_ref_point:
    input:
        "{exp}/profile/{genome}/{bed}/{ref, (TSS)|(TES)|(center)}/{exp}.{samples}.{bed}.{ref}.matrix.gz"
    output:
        "{exp}/profile/{genome}/{bed}/{ref, (TSS)|(TES)|(center)}/{exp}.{samples}.{bed}.{ref}.profile.png"
    params:
        kmeans = "", #"--kmeans 8",
        plotHeight = "5",
        plotWidth = "8"
    threads: 1
    shell:
        "plotProfile --matrixFile {input} --outFileName {output}"
        "{params.kmeans} "
        "--plotHeight {params.plotHeight} "
        "--plotWidth {params.plotWidth} "
        "--refPointLabel {wildcards.ref} "
        "--plotTitle '{wildcards.bed} around {wildcards.ref}: {wildcards.exp}' "

rule profiler_body:
    input:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.matrix.gz"
    output:
        "{exp}/profile/{genome}/{bed}/body/{exp}.{samples}.{bed}.body.profile.png"
    params:
        kmeans = "", #"--kmeans 8",
        sortRegions = "descend", #{descend,ascend,no}
        sortUsing = "mean", #{mean,median,max,min,sum,region_length}
        averageTypeSummaryPlot = "mean", #{mean,median,min,max,std,sum}
        plotHeight = "5",
        plotWidth = "8"
    threads: 1
    shell:
        "plotProfile --matrixFile {input} --outFileName {output} "
        "{params.kmeans} "
        "--plotHeight {params.plotHeight} "
        "--plotWidth {params.plotWidth} "
        "--plotTitle 'body of {wildcards.bed} on {wildcards.exp}' "


########################### NOT WORKING ################################

###############################
#     GC Bias correction      #
###############################
rule computeGCBias_narrow:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        "{exp}/peak_calling/{genome}/{samples}_peaks.narrowPeak",
        "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output:
        "{exp}/mapping/{genome}/GCBias.{samples}.png",
        "{exp}/mapping/{genome}/GCBias.{samples}.matrix.txt"
    #~ threads: 8
    params:
        effectiveGenomeSize = "2451960000",
        genome2bit = lambda wildcards: config["genome"][wildcards.genome]["2bit"]
    log:
        "{exp}/mapping/{genome}/GCBias.{samples}.log"
    run:
        fragmentLength = 200 # default value 
        
        if wildcards.exp in DATA["IP"]:
            fragmentLength = DATA["IP"][wildcards.exp]["fragment_length"]
        
        shell(
            " computeGCBias -b {input[0]} "
            "--genome {params.genome2bit} "
            "--effectiveGenomeSize {params.effectiveGenomeSize} "
            "--filterOut {input[1]} "
            "--fragmentLength " + fragmentLength + " "
            "--numberOfProcessors {threads} "
            "--GCbiasFrequenciesFile {output[1]} "
            "--biasPlot  {output[0]} "
            "2>{log} "
        )

rule computeGCBias_broad:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        "{exp}/peak_calling/{genome}/{samples}_peaks.broadPeak",
        "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output:
        "{exp}/mapping/{genome}/GCBias.{samples}.png",
        "{exp}/mapping/{genome}/GCBias.{samples}.matrix.txt"
    #~ threads: 8
    params:
        effectiveGenomeSize = "2451960000",
        genome2bit = lambda wildcards: config["genome"][wildcards.genome]["2bit"]
    log:
        "{exp}/mapping/{genome}/GCBias.{samples}.log"
    run:
        fragmentLength = 200 # default value 

        if wildcards.exp in DATA["IP"]:
            fragmentLength = DATA["IP"][wildcards.exp]["fragment_length"]

        shell(
            " computeGCBias -b {input[0]} "
            "--genome {params.genome2bit} "
            "--effectiveGenomeSize {params.effectiveGenomeSize} "
            "--filterOut {input[1]} "
            "--fragmentLength " + fragmentLength + " "
            "--numberOfProcessors {threads} "
            "--GCbiasFrequenciesFile {output[1]} "
            "--biasPlot  {output[0]} "
            "2>{log} "
        )

rule computeInputGCBias:
    input:
        "{exp}/mapping/{genome}/{samples}.sorted.bam",
        "{exp}/mapping/{genome}/{samples}.sorted.bam.bai"
    output:
        "{exp}/mapping/{genome}/GCBias.Input.{samples}.png",
        "{exp}/mapping/{genome}/GCBias.Input.{samples}.matrix.txt"
    #~ threads: 8
    params:
        effectiveGenomeSize = "2451960000",
        genome2bit = lambda wildcards: config["genome"][wildcards.genome]["2bit"]
    log:
        "{exp}/mapping/{genome}/GCBias.Input.{samples}.log"
    run:
        fragmentLength = 200 #default
        
        if wildcards.exp in DATA["input"]:
            fragmentLength = DATA["input"][wildcards.exp]["fragment_length"]
        
        shell(
            " computeGCBias -b {input[0]} "
            "--genome {params.genome2bit} "
            "--effectiveGenomeSize {params.effectiveGenomeSize} "
            "--fragmentLength " +fragmentLength + " "
            "--numberOfProcessors {threads} "
            "--GCbiasFrequenciesFile {output[1]} "
            "--biasPlot  {output[0]} "
            "2>{log} "
        )


