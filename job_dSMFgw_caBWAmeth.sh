#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o dSMFalign-output__%I.txt
#BSUB -e dSMFalign-error__%I.txt
#BSUB -J dSMF_caBWAmeth
#BSUB -u jennifer.semple@izb.unibe.ch
#BSUB -N
##BSUB -n 4
##BSUB -R "span[ptile=4]"
#BSUB –R "rusage[mem=8192]" ## in Mb
#BSUB -M 8388608 ## in kb
#BSUB -J array[1]


module add UHTS/Quality_control/fastqc/0.11.5      #fastqc
module add UHTS/Quality_control/cutadapt/1.13     #cutadapt
module add UHTS/Analysis/trimmomatic/0.36; 	#trimmomatic
module add UHTS/Aligner/bwa/0.7.15                 #bwa
#need to install bwa-meth in home directory
module add UHTS/Analysis/samtools/1.4             #samtools
module add UHTS/Analysis/picard-tools/2.9.0        #picard.jar
module add UHTS/Quality_control/qualimap/2.2.1    #qualimap.jar
module add R/3.4.2; #Rscript
module add UHTS/Analysis/BEDTools/2.26.0;	#bedtools	

#before running script, make directory on scratch
runDir=/scratch/cluster/monthly/jsemple/20171214_dSMFgw_N2v182_BWAmeth_cutadapt
scriptDir=/home/jsemple/20171214_dSMFgw_N2v182_BWAmeth_cutadapt
dataDir=/home/jsemple/archive/20171214_dSMFgw_N2v182
#go to the scratch directory and copy script from home directory and run from there
#mkdir -p ${runDir}
#cp -r ${scriptDir}/* ${runDir}/
#cd ${runDir}
#mkdir -p ${runDir}/rawData/
#cp -r ${dataDir}/* ${runDir}/rawData/

bnames=( 180126_SNK268_A_L001_JIB-{1..4} )
let i=${LSB_JOBINDEX}-1
thisbname=${bnames[$i]}

make  bname=${thisbname} all 
