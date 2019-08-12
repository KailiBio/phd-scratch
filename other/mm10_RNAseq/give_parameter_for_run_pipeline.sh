#!/bin/bash

# -- Kaili
# This script is for downloading fastq file and run RNAseq pipeline.

num=$1

#######
workDir="/data/zusers/fankaili/ccre/mm10_rnaseq/"
cd ${workDir}

#
expID=`awk -v num="$num" '{if(NR==num){print $1}}' mouse_RNA_fastq_list_merged.txt`
fileIDs=`awk -v num="$num" '{if(NR==num){print $2}}' mouse_RNA_fastq_list_merged.txt`

# make folder, set path
if [ ! -d ${workDir}rna_exp_Jul25/${expID} ]; then mkdir ${workDir}rna_exp_Jul25/${expID}; fi
cd ${workDir}rna_exp_Jul25/${expID}

if [ ! -d ${workDir}rna_exp_Jul25/${expID}/${fileIDs} ];then mkdir ${workDir}rna_exp_Jul25/${expID}/${fileIDs}; fi

# run pipeline
STARgenomeDir="/home/fankaili/STARgenome"
RSEMrefDir="/home/fankaili/RSEMgenome/RSEMref_mm10"
dataType="unstr_SE"
#
echo bash /data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/run_RNAseq_pipeline.sh ${workDir}rna_exp_Jul25/${expID}/${fileIDs}.fastq "" ${STARgenomeDir} ${RSEMrefDir} ${dataType} 8 8 ${workDir}rna_exp/${expID}/${fileIDs}/ ${fileIDs}
bash /data/zusers/fankaili/github/weng-lab/Kaili/other/mm10_RNAseq/run_RNAseq_pipeline.sh /data/zusers/fankaili/ccre/mm10_rnaseq/fastq/${fileIDs}.fastq "" ${STARgenomeDir} ${RSEMrefDir} ${dataType} 8 8 ${workDir}rna_exp_Jul25/${expID}/${fileIDs}/ ${fileIDs}


# ----
echo "Done. Cheers!"
