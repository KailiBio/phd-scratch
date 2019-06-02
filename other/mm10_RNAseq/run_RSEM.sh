#!/bin/bash

# -- Kaili
# This script is for running RSEM to calculate gene expression.

num=$1

SECONDS=0
#######
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"
cd ${workDir}

#
expID=`awk -v num="$num" '{if(NR==num){print $1}}' mouse_RNA_fastq_list_merged.txt`
fileIDs=`awk -v num="$num" '{if(NR==num){print $2}}' mouse_RNA_fastq_list_merged.txt`

cd ${workDir}rna_exp/${expID}/${fileIDs}/

RSEMrefDir="/home/fankaili/RSEMgenome/RSEMref_mm10"
dataType="unstr_SE"
nThreadsSTAR=8
nThreadsRSEM=8
prefix=${fileIDs}

# executables
RSEM=rsem-calculate-expression

######### RSEM
#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic

# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac

###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir ${prefix} >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir ${prefix} >& Log.rsem

#--------------
echo "5"

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model ${prefix} ${prefix}.pdf
rsem-plot-model ${prefix} ${prefix}.pdf

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

# ----
echo "Done. Cheers!"
