
'''

:date: October 14, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Use HTSEQ-count to count read under exon N-1  / intron / exon / intron /exon N+1

'''
import argparse,textwrap
import subprocess
import logging
from logging.handlers import RotatingFileHandler
from utility import custom_parser
import pysam
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
from datetime import datetime
import os
import HTSeq
import itertools
import numpy
import matplotlib.pylab as plt
from deeptools import getFragmentAndReadSize 

import re
###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################
def write_subprocess_log(completedProcess,logger):
    """
    Write in log the stdout or stderr of subprocess.
    Tcheck if everything was ok.
  
    Args:
        completedProcess (obj): Instance of CompletedProcess send by subprocess.run().
        logger (obj): Instance of logging().
  
    """
    try :
        completedProcess.check_returncode()
        logger.info(completedProcess.stdout)
    except subprocess.CalledProcessError as exc:
                logger.error("===> Exception Caught : ")
                logger.error(exc)  
                logger.error("====> Standard Error : ")
                logger.error(completedProcess.stderr) 
                

def create_logger(config,LEVEL):
    """
    Define a logger instance to write in a file and stream.
  
    Args:
        config (obj): Configuration instance.

    Returns:
        logger (obj): Logger instance to log messages in file and output mainstream.
    
        
    """
    
    logger = logging.getLogger()
    logger.setLevel(LEVEL)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+'activity.log', 'a', 1000000, 1)
    file_handler.setLevel(LEVEL)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(LEVEL)
    logger.addHandler(stream_handler)

    return logger


  
#python3 reshapeBedGraph.py --config=/home/jean-philippe.villemin/code/RNA-SEQ/configs/GHRC38/Heatmap/TEST.json --bed=/home/jean-philippe.villemin/data/data/test.bed###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

    
if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    This script will count signal from bam files. Bam files path are written in json configuration file.  
    Example : 
    python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/htseq_countsignal.py --config=/home/jean-philippe.villemin/code/configs/HEATMAP_HISTONE/K27AC.json

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')

    parameters = parser.parse_args()
    
    config = custom_parser.Configuration(parameters.file_config,"json")

    logger = create_logger(config,logging.INFO)

    logger.info("winInExon "+config.parameters["winInExon"])
    logger.info("winInIntron "+config.parameters["winInIntron"])
    logger.info("list_bed "+",".join(config.parameters["list_bed"].keys()))
    logger.info("list_profile "+",".join(config.parameters["list_profile"]))
    logger.info("path_to_output "+config.parameters["path_to_output"])
    logger.info("files "+",".join(config.parameters["files"].keys()))
    logger.info("plotIs "+config.parameters["plotIs"])

    
    winInExon     = int(config.parameters["winInExon"])
    winInIntron   = int(config.parameters["winInIntron"])
    list_bed      = config.parameters["list_bed"].keys()
    list_profile  = config.parameters["list_profile"]
    
    profiles     = {}
  
    for file2Treat in config.parameters["files"].keys()  :
        
        logger.info(str(file2Treat))
        #if(file2Treat != "K27AC_T7" ) : continue
        profiles[file2Treat]     = {}

        for oneProfile in list_profile : 
            
            profiles[file2Treat][oneProfile] = {}
            profiles[file2Treat][oneProfile]["fragmentSize"]       = int(config.parameters["files"][file2Treat][oneProfile]["fragmentSize"])

            if(config.parameters["duplicates"]=="YES") :
                logger.info("Duplicates is set yes")
                profiles[file2Treat][oneProfile]["bam"]                    =  config.parameters["files"][file2Treat][oneProfile]["bam"] 
                profiles[file2Treat][oneProfile]["TotalMappedReads"]       = pysam.AlignmentFile(config.parameters["files"][file2Treat][oneProfile]["bam"]).mapped
 
            if(config.parameters["duplicates"]=="NO") :
                logger.info("Duplicates is set no")
                profiles[file2Treat][oneProfile]["bam"]                    =  config.parameters["files"][file2Treat][oneProfile]["bam-DupRemoved"] 
                profiles[file2Treat][oneProfile]["TotalMappedReads"]       = pysam.AlignmentFile(config.parameters["files"][file2Treat][oneProfile]["bam-DupRemoved"]).mapped
                
        joinedPathToBed=[]
        for thebedIn in sorted(config.parameters["list_bed"].keys()) :
            
            joinedPathToBed.append(config.parameters["list_bed"][thebedIn]["path"])
            
        godeepTools = "/home/jean-philippe.villemin/code/RNA-SEQ/bash/profileDeeptools.sh "+str(profiles[file2Treat]["Rep1"]["bam"])+" "+str(profiles[file2Treat]["Rep2"]["bam"])+" "+str(profiles[file2Treat]["Control"]["bam"])+" "+str(profiles[file2Treat]["Rep1"]["TotalMappedReads"])+" "+str(profiles[file2Treat]["Rep2"]["TotalMappedReads"])+" "+str(profiles[file2Treat]["Control"]["TotalMappedReads"])+" "+config.parameters["path_to_output"]+" "+",".join(joinedPathToBed)+" "+str(winInExon)+" "+str(winInIntron)+" "+str(file2Treat)+" "+ str(profiles[file2Treat]["Rep1"]["fragmentSize"]) +" "+config.parameters["duplicates"]
        logger.info(godeepTools)
        godeepToolsExec = subprocess.run((godeepTools),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(godeepToolsExec,logger)
        #quit()

        
    
